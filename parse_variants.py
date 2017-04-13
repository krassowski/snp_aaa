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


OFFSET = 20


def show_pos_with_context(seq, start, end):
    return seq[:start] + '→' + seq[start:end] + '←' + seq[end:]


class UnknownChromosome(Exception):
    pass


def get_reference_seq(variant, dna_db, transcript_databases, offset):

    chrom = variant.chr_name

    try:
        chrom = dna_db[chrom]
        seq = chrom.fetch(variant.chrom_start, variant.chrom_end, offset)
    except KeyError:
        raise UnknownChromosome(chrom)

    reference_sequences = []
    for transcript in variant.affected_transcripts:
        ref = get_reference_by_transcript(transcript, transcript_databases, offset)
        reference_sequences.append(ref)

    if not reference_sequences:
        return False

    first = reference_sequences[0]
    if not all(first == ref for ref in reference_sequences):
        print('Potentially different references for different transcripts for %s' % variant.refsnp_id)

    reference_seq = first
    reference_seq['chrom'] = seq

    return reference_seq


def get_reference_by_transcript(transcript, databases, offset):
    reference_seq = {}

    transcript_id = transcript.ensembl_id
    strand = int(transcript.strand)

    for db in databases:
        src = db.sequence_type

        try:
            start = getattr(transcript, src + '_start')
            end = getattr(transcript, src + '_end')
        except IndexError:
            continue

        seq = db.fetch(transcript_id, strand, start, end, offset)

        reference_seq[src] = seq

    return reference_seq


def check_sequence_consistence(reference_seq, ref, alt, offset=OFFSET):
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


def find_sequence(variant, cds_db, cdna_db, dna_db, offset):
    try:
        reference_sequences = get_reference_seq(
            variant,
            dna_db,
            (cdna_db, cds_db),
            offset
        )
    except UnknownChromosome as e:
        print('Unknown chromosome: ', e.message)
        return

    if not reference_sequences:
        return False

    best_sequence = choose_best_seq(reference_sequences)

    if not best_sequence:
        return

    return best_sequence


def determine_mutation(variant, vcf_parser, offset):
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
                    'VCF data contains more than one record matching %s '
                    'variant: %s ' % (variant.refsnp_id, vcf_data)
                )
            elif len(vcf_data) == 0:
                raise ParsingError('No VCF data for %s.' % variant.refsnp_id)

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

    # check ref sequence
    if variant.sequence:
        seq_ref = variant.sequence[offset:-offset]

        if ref != seq_ref:
            print(
                '%s says ref is %s, but sequence analysis pointed to %s for %s'
                % (alt_source, ref, seq_ref, variant.refsnp_id)
            )

    return {
        'gene': gene,
        'ref': ref,
        'chrom_start': pos,
        'alts': alts
    }


def analyze_variant(variant, vcf_parser, cds_db, cdna_db, dna_db, offset=OFFSET):
    variant.correct = True  # at least now

    variant.chrom_start = int(variant.chrom_start)
    variant.chrom_end = int(variant.chrom_end)

    if not variant.sequence:
        surrounding_sequence = find_sequence(variant, cds_db, cdna_db, dna_db, offset)
        if surrounding_sequence:
            variant.sequence = surrounding_sequence
        else:
            print('Cannot determine surrounding sequence for %s variant' % variant)
            # variant.complete = False

    # print('Context: ' + show_pos_with_context(seq, offset, -offset))

    if not (variant.ref and variant.alts and variant.gene and variant.chrom_start):

        mutation_data = determine_mutation(variant, vcf_parser, offset)

        if not mutation_data:
            variant.correct = False
            return

        for attr, value in mutation_data.iteritems():
            old_value = getattr(variant, attr)

            if old_value and old_value != value:
                print(
                    'Deduced %s: %s differs from nominal (%s) for %s'
                    %
                    (attr, value, old_value, variant)
                )

            setattr(variant, attr, value)

    return True


def analyze_poly_a(variant, offset=OFFSET):

    ref_seq = variant.sequence

    # print(variant.refsnp_id)
    # print('Referen: ' + show_pos_with_context(ref_seq, offset, -offset))
    # print('Mutated: ' + show_pos_with_context(mutated_seq, offset, -offset))

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

    return variant


def create_databases():
    constructors = [create_cds_db, create_cdna_db, create_dna_db]

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

    create_dna_db.load_or_create()

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
        '2. Rejection of %s invalid variants' %
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
    unique_variants = set()

    for variant in variants:
        if variant not in unique_variants:
            unique_variants.add(variant)
        else:

            for old in unique_variants:
                if old == variant:
                    old.affected_transcripts.update(variant.affected_transcripts)
                    break

    return unique_variants

