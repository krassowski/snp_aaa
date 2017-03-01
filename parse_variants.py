# -*- coding: utf-8 -*-
from __future__ import print_function
import gc
import re
import signal
import sys
import traceback
from collections import defaultdict
from multiprocessing import Pool

from tqdm import tqdm
from poly_a import poly_a
from snp_parser import vcf_locations
from variant import PolyAAAData
from vcf_parser import ParsingError, VariantCallFormatParser
from snp_parser import create_dna_db, create_cdna_db, create_cds_db


def show_pos_with_context(seq, start, end):
    return seq[:start] + '→' + seq[start:end] + '←' + seq[end:]


class UnknownChromosome(Exception):
    pass


def get_reference_seq(variant, databases, offset):

    reference_seq = {}

    transcript_id = variant.ensembl_transcript_stable_id
    strand = int(variant.ensembl_transcript_chrom_strand)

    for db in databases:

        if type(db) is dict:
            src = 'chrom'
        else:
            src = db.sequence_type

        try:
            start = getattr(variant, src + '_start')
            end = getattr(variant, src + '_end')
        except IndexError:
            #  print(
            #     'Lack of', src, 'coordinates for variant:',
            #     variant.refsnp_id, 'in context of',
            #     transcript_id
            # )
            continue

        if src == 'chrom':
            try:
                chrom = db[variant.chr_name]
                seq = chrom.fetch(start, end, offset)
            except KeyError:
                raise UnknownChromosome(variant.chr_name)
        else:
            seq = db.fetch(transcript_id, strand, start, end, offset)

        reference_seq[src] = seq

    return reference_seq


def check_sequence_consistence(reference_seq, ref, alt, offset=20):
    """Check if CDS and CDNA sequences are consistent with genomic sequence

    (after removal of dashes from both sides).
    """
    ref_seq = reference_seq[ref]
    alt_seq = reference_seq[alt]

    if not alt_seq or ref_seq:
        return

    seq_len = len(alt_seq)
    if not alt_seq or not ref_seq:
        print(
            'Lack of one of required sequences to '
            'check consistence: %s or %s' % (ref, alt)
        )
        return False

    # remove all '-' signs from the beginning of the sequence
    alt_seq = alt_seq.lstrip('-')
    # remember how much was trimmed on the left
    trim_left = seq_len - len(alt_seq)

    seq_len = len(alt_seq)
    alt_seq = alt_seq.rstrip('-')
    trim_right = len(alt_seq) - seq_len

    ref_seq = ref_seq[trim_left:trim_right]

    if alt_seq and ref_seq and alt_seq != ref_seq:
        print('Sequence "%s" is not consistent with sequence "%s":' % (ref, alt))
        print(ref + ':')
        if ref_seq:
            print(
                show_pos_with_context(ref_seq, offset, -offset)
            )
        print(alt + ':')
        if alt_seq:
            print(
                show_pos_with_context(alt_seq, offset, -offset)
            )
        print('')


def choose_best_seq(reference_seq):

    # step 0: check if CDS and CDNA sequences are consistent with genomic seq.
    check_sequence_consistence(reference_seq, 'chrom', 'cds')
    check_sequence_consistence(reference_seq, 'chrom', 'cdna')

    # It's not a big problem if CDS and CDNA sequences are not identical.
    # It's expected that the CDNA will be longer but CDS will be
    # preferred generally.
    # The problem is if they are completely different!
    check_sequence_consistence(reference_seq, 'cds', 'cdna')

    if reference_seq['cdna']:
        chosen = 'cdna'
    elif reference_seq['cds']:
        chosen = 'cds'
    elif reference_seq['chrom']:
        chosen = 'chrom'
    else:
        return False

    return reference_seq[chosen]


def decode_hgvs_code(code):
    """
    Comply to HGVS recommendations: http://www.hgvs.org/mutnomen/recs.html

    tests = [
        'FANCD1:c.658_659delGT',
        'FANCD1:c.5682C>G',
        'FANCD1:c.6275_6276delTT',
        'FANCD1:c.8504C>A',
        'FANCD1:c.5609_5610delTCinsAG',
        'FANCD1:c.1813dupA',
        'FANCD1:c.8219T>A',
        'FANCD1:c.9672dupA'
    ]

    Returns:
        gene, pos, ref, alt
    """
    # TODO testy
    ref, alt = '', ''
    gene, location = code.split(':')
    pos_type = location[0]

    match = re.match(
        '([\d]+)(_[\d]+)?([ACTG]+)?(dup|>|del|ins)([ACTG]+)(dup|>|del|ins)?([ACTG]+)?',
        location[2:]
    )

    if not match:
        raise ValueError('Cannot understand mutation code: %s' % code)

    # genomic and mitochondrial positions can be validated easily
    if pos_type in ('g', 'm'):
        pos = match.group(1)
    elif pos_type in ('c', 'n', 'p'):
        # TODO
        pos = int(match.group(1))
    else:
        raise ParsingError(
            'Wrong type of variant position specification: %s' % pos_type
        )

    event = match.group(4)
    second_event = match.group(6)

    if second_event:
        if event == 'del' and second_event == 'ins':
            ref = match.group(5)
            alt = match.group(7)
        else:
            raise ParsingError(
                'Unknown complicated event: %s and then %s' % (
                    event,
                    second_event
                )
            )
    else:
        if event == '>':
            ref = match.group(3)
            alt = match.group(5)
        elif event in ('dup', 'ins'):
            alt = match.group(5)
        elif event == 'del':
            ref = match.group(5)

    return gene, pos, ref, alt


def analyze_variant(variant, vcf_parser, cds_db, cdna_db, dna_db, offset=20):
    variant.correct = True  # at least now

    variant.chrom_start = int(variant.chrom_start)
    variant.chrom_end = int(variant.chrom_end)

    #print('Checking out variant: %s from %s at %s' % (variant.refsnp_id, variant.refsnp_source, pos))

    try:
        reference_sequences = get_reference_seq(
            variant,
            (cdna_db, cds_db, dna_db),
            offset
        )
    except UnknownChromosome as e:
        print('Unknown chromosome: ', e.message)
        variant.correct = False
        return

    variant.sequence = choose_best_seq(reference_sequences)

    if not variant.sequence:
        variant.correct = False
        print('Skipping: No sequence for ', variant.refsnp_id, 'variant.')

    seq_ref = variant.sequence[offset:-offset]
    #print('Context: ' + show_pos_with_context(seq, offset, -offset))

    if variant.refsnp_source == 'PhenCode':
        alt_source = 'PhenCode'
        gene, pos, ref, alt = decode_hgvs_code(variant.refsnp_id)
        alts = [alt]
    else:
        alt_source = 'VCF'

        try:
            vcf_data = vcf_parser.get_by_variant(variant)
            if len(vcf_data) > 1:
                print(
                    'VCF data containts more than one record matching %s'
                    'variant: %s ' % (variant.refsnp_id, vcf_data)
                )
            analysed_vcf_record = vcf_data[0]
            pos, ref, alts = vcf_parser.parse(analysed_vcf_record)
            gene = vcf_parser.get_gene(analysed_vcf_record)
        except ParsingError as e:
            print(
                'Skipping variant: %s from %s:' % (
                    variant.refsnp_id,
                    variant.refsnp_source
                ),
                end=' '
            )
            print(e.message)
            variant.correct = False
            return False

    variant.alts = alts
    variant.ref = seq_ref

    # HGNC Gene:
    # this causes a problem with CNV mappings
    # should I use that at all?
    if gene:
        variant.gene = gene
    else:
        print('Neither GENEINFO nor GENE are available for', variant.refsnp_id)
        variant.correct = False
        return

    # check ref sequence
    if ref != seq_ref:
        print(
            '%s says ref is %s, but sequence analysis pointed to %s for %s'
            % (alt_source, ref, seq_ref, variant.refsnp_id)
        )

    # check pos
    if variant.chrom_start != pos:
        print(
            'Positions are not matching between %s file and biomart for %s '
            'with ref: %s and alt: %s where positions are: biomart: %s vcf: %s'
            % (
                alt_source, variant.refsnp_id, ref, alts,
                variant.chrom_start, pos
            )
        )

    analyze_poly_a(variant, offset)

    return True


def analyze_poly_a(variant, offset):

    ref_seq = variant.sequence

    #print(variant.refsnp_id)
    #print('Referen: ' + show_pos_with_context(ref_seq, offset, -offset))
    #print('Mutated: ' + show_pos_with_context(mutated_seq, offset, -offset))

    has_aaa, before_len = poly_a(
        ref_seq,
        offset,
        len(ref_seq) - offset
    )

    variant.poly_aaa = defaultdict(PolyAAAData)

    for alt in variant.alts:

        mutated_seq = ref_seq[:offset] + str(alt) + ref_seq[-offset:]

        will_have, after_len = poly_a(
            mutated_seq,
            offset,
            len(mutated_seq) - offset
        )

        variant.poly_aaa[alt].has = has_aaa
        variant.poly_aaa[alt].will_have = will_have
        variant.poly_aaa[alt].before = before_len
        variant.poly_aaa[alt].after = after_len


def create_databases():
    constructors = [create_dna_db, create_cds_db, create_cdna_db]

    databases = [
        db.load()
        for db in constructors
        ]

    return databases


def parse_gene_variants(item):
    gene, variants = item

    # Just to be certain
    variants_unique_ids = set(variant.refsnp_id for variant in variants)
    if len(variants_unique_ids) != len(variants):
        raise Exception('Less unique ids than variants!')

    print('Analysing:', len(variants), 'from', gene)

    vcf_parser = VariantCallFormatParser(
        vcf_locations,
        default_source='ensembl'
    )

    databases = create_databases()

    def analyze_variant_here(variant):

        try:
            analyze_variant(
                variant,
                vcf_parser,
                *databases
            )
        except Exception:
            traceback.print_exc()
            variant.correct = False

        #del databases

        return variant

    variants = map(analyze_variant_here, variants)
    correct_variants = filter(lambda variant: variant.correct, variants)

    print('Analysed', gene)

    sys.stdout.flush()

    gc.collect()

    return gene, correct_variants


def init_worker():
    """Initialization with this function allows to terminate app

    with Ctrl+C even when multiprocessing computation is ongoing.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def parse_variants_by_gene(variants_by_gene):
    """Parses variants"""
    variants_count = sum(
        len(variants) for variants in
        variants_by_gene.values()
    )

    print('Parsing %s variants' % variants_count)

    parsed_variants_by_gene = defaultdict(list)
    merged_duplicates = 0
    rejected_count = 0

    # parsing variants start
    parsing_pool = Pool(12, init_worker, maxtasksperchild=1)

    for gene, variants in tqdm(parsing_pool.imap_unordered(
            parse_gene_variants,
            variants_by_gene.items()
    ), total=len(variants_by_gene)):

        variants_before_parsing = len(variants_by_gene[gene])
        non_unique_parsed = len(variants)

        variants = get_unique_variants(variants)
        unique_parsed = len(variants)

        parsed_variants_by_gene[gene].extend(variants)

        merged_duplicates += non_unique_parsed - unique_parsed
        rejected_count += variants_before_parsing - non_unique_parsed

    parsing_pool.close()
    # parsing end

    parsed_variants_count = sum(
        len(variants) for variants in
        parsed_variants_by_gene.values()
    )

    print(
        '%s variants remained after parsing - quantity was reduced by:\n'
        '1. Merging %s non-unique variants (i.e. identical variants '
        'annotated for different transcripts; you can still find those '
        'transcripts in \'affected_transcripts\' property)\n'
        '2. Rejection of %s variants' %
        (parsed_variants_count, merged_duplicates, rejected_count)
    )

    return parsed_variants_by_gene


def get_unique_variants(variants):
    """Get only unique variants, with respect to:
        - position in genome,
        - set of alternative alleles
        - ensembl gene id

    All transcript affected by variants sharing listed properties
    will expand "affected_transcripts" set of returned variant.
    """
    unique_variants = {}

    for variant in variants:
        transcript = variant.ensembl_transcript_stable_id

        key = '\t'.join(map(
            str,
            [
                variant.chr_name,
                variant.chrom_start,
                variant.chrom_end,
                variant.ref,
                ','.join([str(n) for n in sorted(set(variant.alts))]),
                variant.ensembl_gene_stable_id
            ]
        ))

        if key not in unique_variants.keys():
            unique_variants[key] = variant

        unique_variants[key].affected_transcripts.add(transcript)

    return unique_variants.values()

