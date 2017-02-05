#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from collections import defaultdict
import os
import sys
import gc
import re
import traceback
import cPickle as pickle
import signal
from multiprocessing import Pool
from tqdm import tqdm
from tqdm import trange
from fasta_sequence_db import SequenceDB
from fasta_sequence_db import FastSequenceDB
from cache import cached
from output_formatter import OutputFormatter
from biomart_data import BiomartData
from biomart_data import BiomartDataset
from cna_by_transcript import CompleteCNA
from berkley_hash_set import BerkleyHashSet
from poly_a import poly_a
from variant import Variant
from variant import PolyAAAData
from vcf_parser import ParsingError
from vcf_parser import VariantCallFormatParser


o = OutputFormatter()


GRCH_VERSION = 'GRCh37'
GRCH_SUBVERSION = '13'
ENSEMBL_VERSION = '75'
COSMIC_VERSION = '79'
DBSNP_VERSION = '149'


class VariantsData(BiomartData):

    def __init__(self, dataset=None, filters=None):

        if filters is None:
            filters = {}

        filters.update({
            'so_parent_name':
                [
                    'protein_altering_variant',   # inframe + frameshift
                    'synonymous_variant',
                    'missense_variant',
                    'frameshift_variant'
                ]
        })
        super(self.__class__, self).__init__(dataset, Variant, filters)


def show_pos_with_context(seq, start, end):
    return seq[:start] + '→' + seq[start:end] + '←' + seq[end:]


def gene_names_from_patacsdb_csv(how_many=False):
    """
    Load gene names from given csv file. The default file was exported from
    sysbio.ibb.waw.pl/patacsdb database, from homo sapiens dataset.

    The gene names are expected to be located at second column in each row.
    The count of gene names to read can be constrained with use of an optional
    `how_many` argument (default False denotes that all gene names
    from the file should be returned)
    """
    genes = set()
    count = 0
    with open('patacsdb_all.csv', 'r') as f:
        for line in f:
            if how_many and count >= how_many:
                break
            count += 1
            if line:
                data = [x.strip() for x in line.split(',')]
                genes.add(data[1])

    return list(genes)


def get_variants_by_genes(dataset, gene_names, step_size=50):
    """Retrive from Ensembl's biomart all variants that affect given genes.

    Variants that occur repeatedly in a given gene (i.e. were annotated for
    multiple transcripts) will be reported only once.
    The genes should be specified as the list of gene names.
    """

    variants_by_gene = defaultdict(list)
    variants_count = 0

    print('Downloading variants data from Enseml\'s Biomart:')
    for start in trange(0, len(gene_names), step_size):

        filters = {'ensembl_gene': gene_names[start:start + step_size]}

        print('Downloading variants for genes:', gene_names[start:start + step_size])
        variants = VariantsData(filters=filters, dataset=dataset)
        print('Download completed. Parsing...')

        variants_in_gene = 0

        for variant in variants:
            assert variant.refsnp_id
            variants_in_gene += 1

            gene = variant.ensembl_gene_stable_id

            variants_with_the_same_id = [
                known_variant
                for known_variant in variants_by_gene[gene]
                if known_variant.refsnp_id == variant.refsnp_id
            ]

            if variants_with_the_same_id:
                for known_variant in variants_with_the_same_id:
                    allowed = (
                        # ensembl_transcript_stable_id can differ among two
                        # almost-identical variants
                        variant.ensembl_transcript_stable_id !=
                        known_variant.ensembl_transcript_stable_id
                    )
                    # TODO: check that all other attributes are equal
                    assert allowed
            else:
                variants_by_gene[gene].append(variant)

        # this is a small trick to turn off unpicklable iterator so it is
        # possible to save the variants object as a cache by pickling
        variants.iterator = None

        variants_count += variants_in_gene
        print('Parsed', variants_in_gene, 'variants.')

    print('Downloaded %s variants' % variants_count)
    return variants_by_gene


def get_hgnc(ensembl_transcript_id):

    hgnc_by_ensembl = BerkleyHashSet('hgnc_by_ensembl.db')

    names = hgnc_by_ensembl[ensembl_transcript_id]

    # convert set to a list for indexing support
    names = [name for name in names]

    if len(names) == 0:
        raise ValueError('No HGNC for transcript: ' + ensembl_transcript_id)
    elif len(names) > 1:
        print('Multiple HGNC identifiers for transcript: ' + ensembl_transcript_id)

    return names[0]


class UnknownChromosome(Exception):
    pass


def get_reference_seq(variant, databases, offset):

    reference_seq = {}

    transcript_id = variant.ensembl_transcript_stable_id
    strand = int(variant.ensembl_transcript_chrom_strand)

    for db in databases:
        seq = ''
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


def real_seq_len(sequence):
    return len(sequence.strip('-'))


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
    # The problem is if they are completelly different!
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


vcf_locations = {
    #'COSMIC': 'cosmic/v' + COSMIC_VERSION + '/CosmicCodingMuts.vcf.gz.bgz',
    'dbSNP': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
    GRCH_SUBVERSION + '/00-All.vcf.gz',
    'ClinVar': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
    GRCH_SUBVERSION + '/00-All.vcf.gz',
    'other': 'ensembl/v' + ENSEMBL_VERSION + '/Homo_sapiens.vcf.gz'
}


def decode_phen_code(code):
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
    """
    # TODO testy
    gene, location = code.split(':')
    type = location[0]

    match = re.match(
        '([\d]+)(_[\d]+)?([ACTG]+)?(dup|>|del|ins)([ACTG]+)(dup|>|del|ins)?([ACTG]+)?',
        code
    )

    # genomic and mitohondrial positions can be validated easily
    if type in ('g', 'm'):
        pos = match.group(1)
    elif type in ('c', 'n', 'p'):
        pos = None
    else:
        raise ParsingError(
            'Wrong type of variant position specification: %s' % type
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
            ref = ''
            alt = match.group(5)
        elif event == 'del':
            ref = match.group(5)
            alt = ''

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
        gene, pos, ref, alt = decode_phen_code(variant.refsnp_id)
        alts = [alt]
    else:
        alt_source = 'VCF'

        try:
            vcf_data = vcf_parser.get_by_variant(variant)
            pos, ref, alts = vcf_parser.parse(vcf_data)
            gene = vcf_parser.get_gene(vcf_data)
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


def parse_gene_variants(item):
    gene, variants = item

    # Just to be certain
    variants_unique_ids = set(variant.refsnp_id for variant in variants)
    if len(variants_unique_ids) != len(variants):
        raise Exception('Less unique ids than variants!')

    print('Analysing:', len(variants), 'from', gene)

    vcf_parser = VariantCallFormatParser(vcf_locations)

    constructors = [cachable_dna_db, cachable_cds_db, cachable_cdna_db]

    databases = [
        db()
        for db in constructors
    ]

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


def parse_variants(variants_by_gene):
    """
    vcf.parser uses 0-based coordinates:
    http://pyvcf.readthedocs.org/en/latest/_modules/vcf/parser.html?highlight=coordinates
    ensembl uses 1-based coordinates:
    http://www.ensembl.org/info/docs/api/core/core_tutorial.html#coordinates
    """
    print('Parsing variants:')

    variants_by_gene_by_transcript = {}
    cosmic_genes_to_load = set()

    parsing_pool = Pool(12, init_worker, maxtasksperchild=1)

    for gene, variants in tqdm(parsing_pool.imap_unordered(
            parse_gene_variants,
            variants_by_gene.items()
    ), total=len(variants_by_gene)):

        cosmic_genes_to_load.update(
            [
                variant.gene
                for variant in variants
            ]
        )

        by_transcript = defaultdict(list)

        for variant in variants:
            # The problem with the ensembl's biomart is that it returns records
            # from cosmic without information about the transcript, so we have
            # often a few identical records with only the refsnp_id different,
            # as for example: COSM3391893, COSM3391894
            # Fortunately the transcript id is encoded inside vcf_data retrived
            # from biomart inside the gene identifer (if it is abset, then we
            # have a canonical transcript, at least it is the best guess), eg.:
            # ANKRD26_ENST00000376070 | ANKRD26_ENST00000436985 | ANKRD26
            gene_transcript_id = variant.gene
            transcript_id = gene_transcript_id

            by_transcript[transcript_id].append(variant)

        variants_by_gene_by_transcript[gene] = by_transcript

    parsing_pool.close()

    all_variants_count = sum(
        len(variants) for variants in
        variants_by_gene.values()
    )
    print('All variants', all_variants_count)

    return variants_by_gene_by_transcript, cosmic_genes_to_load


def report(name, data, comment=None):
    """Generates list-based report quickly.

    File will be placed in 'reports' dir and name will be derived from 'name' of
    the report. The data should be an iterable collection of strings.
    Empty reports will not be created.
    """
    directory = 'reports'

    if not os.path.exists(directory):
        os.makedirs(directory)

    if not data:
        return
    with open(directory + '/' + name.replace(' ', '_') + '.txt', 'w') as f:
        f.write('# ' + name + '\n')
        if comment:
            f.write('#' + comment + '\n')
        f.write('\n'.join(data))
    print('Created report "' + name + '" with', len(data), 'entries')


def get_all_variants(variants_by_transcript, gene, report_duplicated=True):

    unique_variants = {}
    duplicated = set()

    for gene_transcript_id, variants in variants_by_transcript.iteritems():
        for variant in variants:
            key = '\t'.join([
                variant.chr_name,
                str(variant.chrom_start),
                str(variant.chrom_end),
                str(variant.ref),
                ','.join([str(n) for n in variant.alts]),
                variant.ensembl_gene_stable_id,
                gene_transcript_id
            ])
            if key not in unique_variants.keys():
                variant.affected_transcripts = set([gene_transcript_id])
                unique_variants[key] = variant
            else:
                unique_variants[key].affected_transcripts.add(gene_transcript_id)
                stored_variant = unique_variants[key]
                duplicated.update([variant.refsnp_id], [stored_variant.refsnp_id])

    if report_duplicated:
        report('Duplicated records for gene: ' + gene, duplicated)

    return unique_variants.values()


def select_poly_a_related_variants(variants):

    return [
        variant
        for variant in variants
        if any([
            data.has or data.will_have
            for data in variant.poly_aaa.values()
        ])
    ]


def summarize_poly_aaa_variants(variants_by_gene_by_transcript):

    variant_aaa_report = []

    for gene, variants_by_transcript in variants_by_gene_by_transcript.iteritems():

        # treat all variants the same way - just remove duplicates
        variants = get_all_variants(variants_by_transcript, gene)

        poly_a_related_variants = select_poly_a_related_variants(variants)

        variant_aaa_report += ['# ' + gene]
        variant_aaa_report += [
            '\t'.join(map(str, [
                variant.refsnp_id,
                data.increased,
                data.decreased,
                data.change,
                alt
            ]))
            for variant in poly_a_related_variants
            for alt, data in variant.poly_aaa.items()
        ]

    report(
        'poly aaa increase and decrease by variants',
        variant_aaa_report,
        'snp_id\tpoly_aaa_increase\tpoly_aaa_decrease\tpoly_aaa_change\talt'
    )


def gtex_over_api(variants_by_gene_by_transcript):
    import requests
    import sys

    variant_aaa_report = []

    for gene, variants_by_transcript in tqdm(variants_by_gene_by_transcript.iteritems()):

        # treat all variants the same way - just remove duplicates
        variants = get_all_variants(variants_by_transcript, gene)
        poly_a_related_variants = select_poly_a_related_variants(variants)

        """
        print(
            'Checikng %s (out of %s) poly A related variants from %s.' %
            (
                len(poly_a_related_variants),
                len(variants),
                gene
            )
        )
        """

        for variant in tqdm(poly_a_related_variants):

            #from expression_database import TISSUES_LIST as tissues

            server = 'http://rest.ensembl.org'
            #for tissue in tissues:
            ext = '/eqtl/variant_name/homo_sapiens/' + variant.refsnp_id + '?statistic=p-value;content-type=application/json'

            #"rs17438086?statistic=p-value;stable_id=ENSG00000162627;tissue="

            u_a = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/48.0.2564.82 Safari/537.36'

            try:
                r = requests.get(
                    server + ext,
                    headers={
                        'USER-AGENT': u_a,
                        'content-type': 'application/json'
                    }
                )

                if not r.ok:
                    r.raise_for_status()
                    sys.exit()

                decoded = r.json()

                if not decoded['error']:
                    print('Found sth!')
                    print(repr(decoded))

            except:
                pass


def summarize_copy_number_expression(variants_by_gene_by_transcript, cna):

    variants_count = 0
    poly_a_related_variants_count = 0

    cnv_aaa_report = []
    import operator

    no_expression = set()

    for gene, variants_by_transcript in variants_by_gene_by_transcript.iteritems():

        expression = [0, 0, 0]
        for gene_transcript_id in variants_by_transcript.keys():

            try:
                expr = cna.get_by_gene_and_transcript(gene_transcript_id)
                expression = map(operator.add, expression, expr)
            except KeyError:
                # no expression data for this transcript
                no_expression.add(gene_transcript_id)
                continue

        # treat all variants the same way - just remove duplicates
        variants = get_all_variants(variants_by_transcript, gene, report_duplicated=False)

        """
        # for this analysis consider only variants from cosmic
        variants = [
            variant
            for variant in variants
            if variant.refsnp_source == 'COSMIC'
        ]
        """

        variants_count += len(variants)

        print(gene, 'with its', len(variants), 'variants analysed')

        # loss of poly_aaa, decrease in length (-1 per residue)
        decrease = 0
        # gain of poly_aaa: , increase in length (+1 per residue)
        increase = 0

        poly_a_related_variants = select_poly_a_related_variants(variants)
        poly_a_related_variants_count += len(poly_a_related_variants)

        # give scores for length increase
        for variant in poly_a_related_variants:
            if variant.poly_aaa.increased:
                increase += 1
            elif variant.poly_aaa.decreased:
                decrease += 1

        cnv_aaa_report += [
            (gene, expression[0], expression[2], increase, decrease)
        ]

    report('no expression data for some transcripts', no_expression)

    gene_count = len(variants_by_gene_by_transcript)

    print('analyzed genes:', gene_count)
    print('analyzed variants:', variants_count)
    print('poly_a related variants', poly_a_related_variants_count)

    report(
        'poly a and expression table',
        ['\t'.join(map(str, line)) for line in cnv_aaa_report],
        'gene\tcnv+\tcnv-\taaa+\taaa-'
    )


def get_all_used_transcript_ids(variants_by_gene):
    """Return all transcript identificators that occur in the passed variants.

    variants_by_gene is expected to by a dict of variant lists grouped by genes.
    """
    transcripts_ids = set()

    for variants in variants_by_gene.itervalues():
        transcripts_ids.update({
            variant.ensembl_transcript_stable_id
            for variant in variants
        })

    return transcripts_ids


def poly_aaa_vs_expression(variants_by_gene_by_transcript):

    from expression_database import ExpressionDatabase

    bdb = ExpressionDatabase('expression_slope_in_tissues_by_mutation_full.db')

    def is_length_difference_big(l1, l2):
        """Is the first list much longer than the second?"""
        len1 = len(l1)
        len2 = len(l2)
        assert len1 > len2

        if len2 == 0 or len1 / len2 > 10:
            return True

    gtex_report = []
    gtex_report_by_genes = []

    for gene, variants_by_transcript in variants_by_gene_by_transcript.iteritems():

        # treat all variants the same way - just remove duplicates
        variants = get_all_variants(variants_by_transcript, gene, report_duplicated=False)
        poly_a_related_variants = select_poly_a_related_variants(variants)

        print('Analysing %s poly_a related vartiants (total: %s) from %s gene.' % (len(poly_a_related_variants), len(variants), gene))

        for variant in poly_a_related_variants:

            variant.expression = {}
            expression_data_by_alt = bdb.get_by_mutation(variant)

            for alt, expression_data in expression_data_by_alt.items():

                if not expression_data:
                    # print('No expression for', variant.refsnp_id)
                    expression_trend = 'none'
                    continue
                else:
                    print('Expression data for', variant.refsnp_id, 'found:', expression_data)

                expression_up = []
                expression_down = []

                for tissue, slope in expression_data:
                    if slope > 0:
                        expression_up += [tissue]
                    elif slope < 0:
                        expression_down += [tissue]

                # is this rather up?
                if len(expression_up) > len(expression_down):
                    # is this certainly up?
                    if is_length_difference_big(expression_up, expression_down):
                        expression_trend = 'up'
                    else:
                        expression_trend = 'rather_up'
                # is this rather down?
                elif len(expression_down) > len(expression_up):
                    # is this certainly down?
                    if is_length_difference_big(expression_down, expression_up):
                        expression_trend = 'down'
                    else:
                        expression_trend = 'rather_down'
                # is unchanged?
                else:
                    expression_trend = 'constant'

                expression_up_in_X_cases = len(expression_up)
                expression_down_in_X_cases = len(expression_down)

                data = variant.poly_aaa[alt]
                variant.expression[alt] = expression_trend

                gtex_report += [(
                    variant.refsnp_id,
                    expression_up_in_X_cases,
                    expression_down_in_X_cases,
                    expression_trend,
                    data.increased,
                    data.decreased,
                    data.change,
                    alt
                )]

        gtex_report_by_genes += [(
            gene,
            sum('up' in v.expression.values() for v in poly_a_related_variants),
            sum('down' in v.expression.values() for v in poly_a_related_variants),
            sum(
                sum('up' == expr for expr in v.expression.values())
                for v in poly_a_related_variants
            ),
            sum(
                sum('down' == expr for expr in v.expression.values())
                for v in poly_a_related_variants
            ),
            sum(data.increased for v in poly_a_related_variants for data in v.poly_aaa.values()),
            sum(data.decreased for v in poly_a_related_variants for data in v.poly_aaa.values())
        )]

    report(
        'Expression table for variants (based on data from GTEx)',
        ['\t'.join(map(str, line)) for line in gtex_report],
        'variant\texpression+\texpression-\ttrend\taaa+\taaa-\taaa change'
    )

    report(
        'Expression table for genes (based on data from GTEx)',
        ['\t'.join(map(str, line)) for line in gtex_report_by_genes],
        # note: alleles is not the same as variants
        'gene\talleles with expression+\talleles with expression-\tvariants with expression+\tvariants with expression-\t#aaa+\t#aaa-'
    )

    print('Done')


def main(args, dataset):
    """The main workflow happens here."""

    variants_by_gene = cachable_variants_by_gene()

    if args.parse_variants:

        variants_by_gene_by_transcript, cosmic_genes_to_load = parse_variants(
            variants_by_gene
        )

        try:
            with open('.variants_by_gene_by_transcript_37_all_alts.cache', 'wb') as f:
                pickle.dump(variants_by_gene_by_transcript, f, protocol=pickle.HIGHEST_PROTOCOL)
            with open('.cosmic_genes_to_load_37_all_alts.cache', 'wb') as f:
                pickle.dump(cosmic_genes_to_load, f, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception:
            traceback.print_exc()
    else:
        with open('.variants_by_gene_by_transcript_37_all_alts.cache', 'rb') as f:
            variants_by_gene_by_transcript = pickle.load(f)
        with open('.cosmic_genes_to_load_37_all_alts.cache', 'rb') as f:
            cosmic_genes_to_load = pickle.load(f)

    try:
        if 'gtex_over_api' in args.report:
            gtex_over_api(variants_by_gene_by_transcript)
    except Exception:
        traceback.print_exc()

    try:
        if 'list_poly_aaa_variants' in args.report:
            summarize_poly_aaa_variants(variants_by_gene_by_transcript)
    except Exception:
        traceback.print_exc()

    try:
        if 'copy_number_expression' in args.report:
            # TODO to cache again later
            @cached(action='load')
            def cachable_cna():
                return CompleteCNA(
                    'cosmic/v' + COSMIC_VERSION + '/CosmicCompleteCNA.tsv',
                    restrict_to=cosmic_genes_to_load
                )

            cna = cachable_cna()
            summarize_copy_number_expression(variants_by_gene_by_transcript, cna)
    except Exception:
        traceback.print_exc()

    try:
        if 'poly_aaa_vs_expression' in args.report:
            poly_aaa_vs_expression(variants_by_gene_by_transcript)
    except Exception:
        traceback.print_exc()


if __name__ == '__main__':

    run = False

    import argparse

    to_show = ['databases', 'datasets', 'filters', 'attributes',
               'attributes_by_page', 'some_variant']
    parser = argparse.ArgumentParser(description='Find SNPs')
    parser.add_argument('--show', choices=to_show)
    parser.add_argument(
        '--report',
        nargs='+',
        choices=[
            'list_poly_aaa_variants', 'copy_number_expression',
            'poly_aaa_vs_expression', 'gtex_over_api'
        ],
        default=[]
    )
    parser.add_argument('--profile', action='store_true')
    parser.add_argument(
        '--dataset',
        type=str,
        help='name of biomart dataset to be used, eg. hsapiens_snp',
        default='hsapiens_snp_som'
    )
    parser.add_argument(
        '--biomart',
        type=str,
        help='URL of biomart to be used. '
             'For ensembl mirrors replace www with: uswest, useast or asia',
        default='http://www.ensembl.org/biomart'
    )
    parser.add_argument(
        '-n',
        '--number',
        type=int,
        help='Number of genes to analyze',
        default=None
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='increase output verbosity'
    )
    parser.add_argument(
        '-r',
        '--reload_variants',
        action='store_true',
    )
    parser.add_argument(
        '--parse_variants',
        action='store_true',
    )
    parser.add_argument(
        '-g',
        '--reload_gtex',
        action='store_true',
    )
    parser.add_argument(
        '-s',
        '--step_size',
        type=int,
        default=25,
        help='Variants from how many genes should be requested at once from biomart? Default 25.'
    )

    subcommands = ['show', 'profile']
    arguments = ['--' + a if a in subcommands else a for a in sys.argv[1:]]

    args = parser.parse_args(arguments)

    global BIOMART_URL
    BIOMART_URL = args.biomart

    snp_dataset = BiomartDataset(args.biomart, name=args.dataset)

    o.force = args.verbose

    if args.show:
        what = args.show
        if what == 'databases':
            from biomart import BiomartServer
            BiomartServer(BIOMART_URL).show_databases()
        if what == 'datasets':
            from biomart import BiomartServer
            BiomartServer(BIOMART_URL).show_datasets()
        if what == 'filters':
            snp_dataset.show_filters()
        if what == 'attributes':
            snp_dataset.show_attributes()
        if what == 'attributes_by_page':
            snp_dataset.show_attributes_by_page()
        if what == 'some_variant':
            # it would be very strange if we do not find
            # any variants in first 5 random genes
            genes_from_patacsdb = gene_names_from_patacsdb_csv(how_many=5)
            variants_by_gene = get_variants_by_genes(
                snp_dataset,
                genes_from_patacsdb
            )
            gene, variants = variants_by_gene.popitem()
            print(variants[0])
    else:

        if args.profile:
            import profile
            profile.run('main(args, snp_dataset)')
        else:
            run = True

    if args.reload_variants:
        global_cache_action = 'save'
    elif run:
        global_cache_action = 'load'

    if run or args.reload_variants:

        genes_from_patacsdb = gene_names_from_patacsdb_csv(args.number)

        @cached(action='load')
        def cachable_variants_by_gene():
            return get_variants_by_genes(snp_dataset, genes_from_patacsdb, step_size=args.step_size)

        variants_by_gene = cachable_variants_by_gene()


        @cached(action='load')
        def cachable_transcripts_to_load():
            return get_all_used_transcript_ids(variants_by_gene)

        transcripts_to_load = cachable_transcripts_to_load()

        @cached(action=global_cache_action)
        def cachable_cds_db():
            return SequenceDB(
                version=ENSEMBL_VERSION,
                assembly=GRCH_VERSION,
                index_by='transcript',
                sequence_type='cds',
                restrict_to=transcripts_to_load
            )

        @cached(action=global_cache_action)
        def cachable_cdna_db():
            return SequenceDB(
                version=ENSEMBL_VERSION,
                assembly=GRCH_VERSION,
                index_by='transcript',
                sequence_type='cdna',
                restrict_to=transcripts_to_load
            )


        @cached(action=global_cache_action)
        def cachable_dna_db():
            chromosomes = map(str, range(1, 23)) + ['X', 'Y', 'MT']
            dna_db = {}

            for chromosome in chromosomes:
                dna_db[chromosome] = FastSequenceDB(
                    version=ENSEMBL_VERSION,
                    assembly=GRCH_VERSION,
                    sequence_type='dna',
                    id_type='chromosome.' + chromosome
                )
            return dna_db

        constructors = [cachable_dna_db, cachable_cds_db, cachable_cdna_db]

        for db in constructors:
            db()

        print('Variants data ' + global_cache_action + 'ed')

    if args.reload_gtex:
        from expression_database import ExpressionDatabase
        from expression_database import import_expression_data

        bdb = ExpressionDatabase('expression_slope_in_tissues_by_mutation_full.db')

        print('Reloading GTEx expression data:')
        import_expression_data(
            bdb,
            path='GTEx_Analysis_v6p_all-associations',
            suffix='_Analysis.v6p.all_snpgene_pairs.txt.gz'
        )

    if run and not args.reload_variants:

        main(args, snp_dataset)
