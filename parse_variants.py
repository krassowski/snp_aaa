# -*- coding: utf-8 -*-
from __future__ import print_function
import gc
import re
import signal
import traceback
from collections import defaultdict
from multiprocessing import Pool

from pyfaidx import Fasta
from tqdm import tqdm
from poly_a import poly_a
from snp_parser import vcf_mutation_sources, jit
from variant import PolyAAAData
from vcf_parser import ParsingError, VariantCallFormatParser
from snp_parser import create_dna_db
from snp_parser import TRANSCRIPT_DB_PATH


REFERENCE_SEQUENCE_TYPE = 'cds'
OFFSET = 20
# speeds up parsing a lot if one is looking only for poly aaa variants
KEEP_ONLY_POLY_A = True


def show_pos_with_context(seq, start, end):
    return seq[:start] + '>' + seq[start:end] + '<' + seq[end:]


def get_dna_sequence(variant, offset, database):
    chromosome = variant.chr_name
    try:
        chrom_db = database[chromosome]
        return chrom_db.fetch(variant.chrom_start, variant.chrom_end, variant.chrom_strand, offset)
    except KeyError:
        print('Unknown chromosome: %s' % chromosome)


transcripts = Fasta(TRANSCRIPT_DB_PATH, key_function=lambda x: x.split('.')[0])


def get_transcripts_sequences(variant, offset):
    reference_sequences = {}

    any_ref = False

    for transcript in variant.affected_transcripts:

        #ref = get_reference_by_transcript(transcript, offset, database)

        raw_start = getattr(transcript, REFERENCE_SEQUENCE_TYPE + '_start')
        raw_end = getattr(transcript, REFERENCE_SEQUENCE_TYPE + '_end')

        if not (raw_start and raw_end):
            # print('No sequence coordinates for %s transcript' % transcript.ensembl_id)
            continue

        if raw_end - raw_start > 2 * offset:
            print('Skipping transcript %s of variant %s: too wide variation' % (transcript.ensembl_id, variant.refsnp_id))
            continue

        if transcript.ensembl_id not in transcripts:
            continue

        whole_seq = str(transcripts[transcript.ensembl_id])

        seq = None
        if whole_seq:
            start = raw_start - 1
            end = raw_end

            cut_from = start - offset
            cut_to = end + offset

            seq = whole_seq[max(cut_from, 0):max(cut_to, 0)]

            # if we are on the edge of sequence
            if cut_from < 0:
                seq = '-' * (-cut_from) + seq
            if cut_to > len(whole_seq):
                seq += '-' * (cut_to - len(whole_seq))

            assert len(seq) == offset * 2 + raw_end - raw_start + 1

        ref = seq
        #assert ref == seq

        if ref:
            reference_sequences[transcript] = ref
            any_ref = True
        elif transcript.ensembl_id.startswith('LRG_'):
            # Yeah, I know that LRG transcripts do not have sequences in ensembl databases
            pass
        else:
            print('No sequence for %s transcript' % transcript)

    if not any_ref:
        print(
            'Variant %s has no sequence for any of its transcripts (using %s)'
            %
            (variant.refsnp_id, REFERENCE_SEQUENCE_TYPE)
        )

    return reference_sequences


def get_reference_by_transcript(transcript, offset, database):
    src = database.sequence_type

    start = getattr(transcript, src + '_start')
    end = getattr(transcript, src + '_end')

    if start is None or end is None:
        return False

    return database.fetch(transcript.ensembl_id, transcript.strand, start, end, offset)


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


def determine_mutation(variant, vcf_parser, offset):
    if variant.refsnp_source == 'PhenCode':
        print('PhenCode variant:')
        print(variant)
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
                print('Only the first record will be parsed!')
            elif len(vcf_data) == 0:
                raise ParsingError('No VCF data for %s.' % variant.refsnp_id)

            analysed_vcf_record = vcf_data[0]
            strand = list(variant.affected_transcripts)[0].strand
            for transcript in variant.affected_transcripts:
                if strand != transcript.strand:
                    print('Stand mismatch for')
                    print(variant)
            pos, ref, alts = vcf_parser.parse(analysed_vcf_record, variant.refsnp_source, strand)
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
    if variant.sequences:

        pos_correct = 0

        for transcript, sequence in variant.sequences.iteritems():
            seq_ref = sequence[offset:-offset]

            if ref != seq_ref:

                # Cosmic represents insertions AND bigger deletions as a range:
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=111380 - deletion (start = end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=1223920 - substitution (start = end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=4699526 - insertion (start != end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=4612790 - big deletion (start != end)

                if variant.refsnp_source == 'COSMIC' and transcript.cds_end != transcript.cds_start:
                    # as all deletions pass ref == seq_ref test, and insertions fall there,
                    # it is needed to assure that those are corrected here. Firstly, make sure
                    # that it is an insertion:
                    assert all(len(alt) > len(ref) for alt in alts)

                    del variant.sequences[transcript]

                    # sequence was too big by 'diff'
                    diff = transcript.cds_end - transcript.cds_start + 1
                    if not pos_correct:
                        pos_correct = diff - 1
                    else:
                        assert pos_correct == diff - 1

                    sequence = sequence[:-diff]
                    transcript.cds_end = transcript.cds_start
                    assert sequence[offset:-offset] == ref

                    variant.sequences[transcript] = sequence

                else:

                    print(
                        '%s says ref is %s, but sequence analysis pointed to %s for %s, %s'
                        % (alt_source, ref, seq_ref, variant.refsnp_id, transcript.ensembl_id)
                    )
        pos -= pos_correct

    return {
        'gene': gene,
        'ref': ref,
        'chrom_start': pos,
        'alts': alts
    }


def analyze_variant(variant, vcf_parser, offset=OFFSET, keep_only_poly_a=False):
    variant.correct = True  # at least now

    #if database.sequence_type == 'dna':
    #    if variant.chrom_end - variant.chrom_start > 2 * offset:
    #        print('Skipping variant %s: too wide variation' % variant.refsnp_id)
    #        variant.correct = False
    #        return
    #    chosen_surrounding_sequence = get_dna_sequence(variant, offset, database)
    #    sequences = {'dna': chosen_surrounding_sequence}
    #else:
    sequences = get_transcripts_sequences(variant, offset)

    if keep_only_poly_a:
        retained_sequences = {}

        for transcript, sequence in sequences.iteritems():

            accepted, score = poly_a(sequence, offset, len(sequence) - offset)
            # assuming that we have only single, point substitutions
            worst_mutated_seq = sequence[:offset] + 'A' + sequence[-offset:]

            will_have, after_len = poly_a(
                worst_mutated_seq,
                offset,
                len(worst_mutated_seq) - offset
            )

            if accepted or will_have:
                print('Variant with poly A: %s, %s' % (variant.refsnp_id, transcript.ensembl_id))
                print(show_pos_with_context(sequence, offset, -offset))
                retained_sequences[transcript] = sequence

        if not retained_sequences:
            # print('Variant %s skipped (no chances for poly a)' % variant.refsnp_id)
            variant.correct = False
            return

        sequences = retained_sequences

    if sequences:
        variant.sequences = sequences
    else:
        print('Cannot determine surrounding sequences for %s variant' % variant.refsnp_id)
        variant.correct = False
        return

    if not (variant.ref and variant.alts and variant.gene and variant.chrom_start):

        mutation_data = determine_mutation(variant, vcf_parser, offset)

        if not mutation_data:
            variant.correct = False
            return

        for attr, value in mutation_data.iteritems():
            old_value = getattr(variant, attr)

            if old_value:
                if old_value != value:
                    print(
                        'Deduced %s: %s differs from nominal (%s) for %s'
                        %
                        (attr, value, old_value, variant)
                    )
                    setattr(variant, attr, value)
                # else old_value == value, no need for setattr
            else:
                setattr(variant, attr, value)

    return True


def analyze_poly_a(variant, offset=OFFSET):

    for transcript, sequence in variant.sequences.iteritems():
        poly_aaa = get_poly_a(sequence, variant.alts, offset)
        transcript.poly_aaa = poly_aaa

    return variant


def get_poly_a(ref_seq, alts, offset):

    has_aaa, before_len = poly_a(
        ref_seq,
        offset,
        len(ref_seq) - offset
    )

    poly_aaa = defaultdict(PolyAAAData)

    for alt in alts:

        mutated_seq = ref_seq[:offset] + str(alt) + ref_seq[-offset:]

        will_have, after_len = poly_a(
            mutated_seq,
            offset,
            len(mutated_seq) - offset
        )

        poly_aaa[alt].has = has_aaa
        poly_aaa[alt].will_have = will_have
        poly_aaa[alt].before = before_len
        poly_aaa[alt].after = after_len

    return poly_aaa


def parse_gene_variants(item):
    gene, variants = item

    # Just to be certain
    variants_unique_ids = set(variant.refsnp_id for variant in variants)
    if len(variants_unique_ids) != len(variants):
        raise Exception('Less unique ids than variants!')

    print('Analysing:', len(variants), 'variants from some group')

    vcf_parser = VariantCallFormatParser(
        vcf_mutation_sources,
        default_source='ensembl'
    )

    #database = create_database()

    def analyze_variant_here(variant):

        try:
            analyze_variant(variant, vcf_parser, keep_only_poly_a=KEEP_ONLY_POLY_A)
        except Exception:
            traceback.print_exc()
            variant.correct = False

        return variant

    variants = map(analyze_variant_here, variants)

    correct_variants = filter(lambda variant: variant.correct, variants)

    print('Analysed:', len(variants), 'variants from some group')

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
        variants_by_gene.itervalues()
    )

    print('Parsing %s variants' % variants_count)

    parsed_variants_by_gene = defaultdict(list)
    merged_duplicates = 0
    rejected_count = 0

    create_dna_db.load_or_create()

    # parsing variants start
    parsing_pool = Pool(3, init_worker, maxtasksperchild=1)

    for gene, variants in tqdm(parsing_pool.imap_unordered(
            parse_gene_variants,
            variants_by_gene.iteritems()
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
        parsed_variants_by_gene.itervalues()
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


#@jit
def get_unique_variants(variants):
    """Get only unique variants, with respect to:
        - position in genome,
        - set of alternative alleles
        - ensembl gene id

    All transcript affected by variants sharing listed properties
    will expand "affected_transcripts" set of returned variant.
    """
    # adding and removing None is workaround for https://github.com/numba/numba/issues/2152
    unique_variants = set()

    for variant in variants:
        if variant not in unique_variants:
            unique_variants.add(variant)
        else:

            for old in unique_variants:
                if old == variant:
                    old.affected_transcripts.update(variant.affected_transcripts)
                    break

    #unique_variants.remove(None)

    return unique_variants

