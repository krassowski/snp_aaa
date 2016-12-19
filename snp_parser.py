#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from collections import defaultdict
import vcf
import os
import sys
from poly_a import poly_a
from fasta_sequence_db import SequenceDB
from fasta_sequence_db import FastSequenceDB
from cache import cached
from output_formatter import OutputFormatter
from biomart_data import BiomartData
from biomart_data import BiomartDataset
from cna_by_transcript import CompleteCNA
from tqdm import tqdm
from berkley_hash_set import BerkleyHashSet
import gc
import traceback


o = OutputFormatter()


class VariantsData(BiomartData):

    def __init__(self, dataset=None, attributes=None, filters=None):

        if attributes is None:
            attributes = []
        if filters is None:
            filters = {}

        attributes += [
            'refsnp_id',
            'refsnp_source',
            'chr_name',
            'chrom_start',
            'chrom_end',
            'allele',  # Variant Alleles
            'allele_1',  # Ancestral allele
            'minor_allele',
            'chrom_strand',
            'cdna_start',
            'cdna_end',
            'ensembl_gene_stable_id',
            'ensembl_transcript_stable_id',
            'ensembl_transcript_chrom_strand',
            'cds_start',
            'cds_end'
            # 'consequence_type_tv',
            # 'consequence_allele_string'
        ]

        filters.update({
            'so_parent_name':
                [
                    'protein_altering_variant',   # inframe + frameshift
                    'synonymous_variant',
                    'missense_variant',
                    'frameshift_variant'
                ]
        })
        super(self.__class__, self).__init__(dataset, attributes, filters)


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


def get_variants_by_genes(dataset, gene_names):
    """Retrive from Ensembl's biomart all variants that affect given genes.

    Variants that occur repeatedly in a given gene (i.e. were annotated for
    multiple transcripts) will be reported only once.
    The genes should be specified as the list of gene names.
    """

    variants_by_gene = defaultdict(list)

    for start in range(0, len(gene_names), 300):

        filters = {'ensembl_gene': gene_names[start:start + 300]}

        variants = VariantsData(filters=filters, dataset=dataset)

        for variant in variants:

            gene = variant.ensembl_gene_stable_id

            assert variant.refsnp_id

            for known_variant in variants_by_gene[gene]:
                if variant.refsnp_id == known_variant.refsnp_id:
                    allowed = (
                        # ensembl_transcript_stable_id can differ among two
                        # almost-identical variants
                        variant.ensembl_transcript_stable_id !=
                        known_variant.ensembl_transcript_stable_id
                    )
                    # TODO: check that all other attributes are equal
                    assert allowed
                    break
            else:
                variants_by_gene[gene].append(variant)

        # this is a small trick to turn off unpicklable iterator so it is
        # possible to save the variants object as a cache by pickling
        variants.iterator = None

    return variants_by_gene


def get_hgnc(variant, hgnc_by_ensembl):
    return hgnc_by_ensembl[variant.ensembl_transcript_stable_id]


def get_vcf_by_variant(pos, variant):

    if variant.refsnp_source == 'COSMIC':
        vcf_reader = vcf.Reader(filename='cosmic/CosmicCodingMuts.vcf.gz.bgz')
        record_id = variant.refsnp_id
    elif variant.refsnp_source == 'dbSNP':
        vcf_reader = vcf.Reader(filename='ncbi/00-All.vcf.gz')
        record_id = variant.refsnp_id
        print('Using record_id: ' + record_id)
    else:
        vcf_reader = vcf.Reader(filename='ensembl/Homo_sapiens_somatic.vcf.gz')
        hgnc_by_ensembl = BerkleyHashSet('hgnc_by_ensembl.db')

        if variant.cds_start is None:
            raise ValueError(
                'No cds_start data for ' + variant.ensembl_transcript_stable_id
            )

        names = get_hgnc(variant, hgnc_by_ensembl)
        names_list = [name for name in names]
        if len(names_list) > 1:
            print('Multiple HGNC identifiers for transcript: ' + variant.ensembl_transcript_stable_id)
        elif len(names_list) == 0:
            raise ValueError(
                'No HGNC for transcript: ' +
                variant.ensembl_transcript_stable_id
            )
        else:
            name = names_list[0]

        allele = variant.allele.replace('/', '>')
        record_id = ''.join([name, ':c.', variant.cds_start, allele])

    return get_vcf_by_id(vcf_reader, pos, record_id)


def get_vcf_by_id(vcf, pos, record_id):

    vcf_data = None

    for record in vcf.fetch(*pos):
        if record.ID == record_id:
            vcf_data = record
            break

    return vcf_data


def get_reference_seq(variant, databases, offset):

    reference_seq = {}

    transcript_id = variant.ensembl_transcript_stable_id
    strand = int(variant.ensembl_transcript_chrom_strand)

    for db in databases:
        seq = ''
        src = db.sequence_type

        try:
            start = variant.get(src + '_start')
            end = variant.get(src + '_end')
        except IndexError:
            #  print(
            #     'Lack of', src, 'coordinates for variant:',
            #     variant.refsnp_id, 'in context of',
            #     transcript_id
            # )
            continue

        if src == 'dna':
            try:
                chrom = db[variant.chr_name]
                seq = chrom.fetch(start, end, offset)
            except KeyError as e:
                print('Unknown chromosome', variant.chr_name)
                raise e
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

    if alt_seq != ref_seq:
        print('Sequence %s is not consistent sequence %s:' % (ref, alt))
        if ref_seq:
            print(
                ref_seq, ':\t', show_pos_with_context(ref_seq, offset, -offset)
            )
        if alt_seq:
            print(
                alt_seq, ':\t', show_pos_with_context(alt_seq, offset, -offset)
            )


def choose_best_seq(reference_seq):

    # step 0: check if CDS and CDNA sequences are consistent with genomic seq.
    check_sequence_consistence(reference_seq, 'dna', 'cds')
    check_sequence_consistence(reference_seq, 'dna', 'cdna')

    # It's not a big problem if CDS and CDNA sequences are not identical.
    # It's expected that the CDNA will be longer but CDS will be
    # preferred generally.
    # The problem is if they are completelly different!
    check_sequence_consistence(reference_seq, 'cds', 'cdna')

    if reference_seq['cdna']:
        chosen = 'cdna'
    elif reference_seq['cds']:
        chosen = 'cds'
    elif reference_seq['dna']:
        chosen = 'dna'
    else:
        return False

    return reference_seq[chosen]


def analyze_variant(variant, cds_db, cdna_db, dna_db, offset=20):
    variant.correct = True  # at least now

    variant.chrom_start = int(variant.chorm_start)
    variant.chrom_end = int(variant.chorm_end)

    pos = [variant.chr_name, variant.chrom_start, variant.chrom_end]

    variant.length = variant.chrom_end - variant.chrom_start + 1  # TODO test this

    try:
        reference_sequences = get_reference_seq(
            variant,
            (cdna_db, cds_db, dna_db),
            offset
        )
    except KeyError:
        variant.correct = False
        return

    # to get to vcf stored data by vcf reader, change coordinates to 0-based
    pos[1] -= 1
    pos[2] -= 1

    # and represent them as range
    # (if you had n:n pointing to a single base, use n:n+1)
    pos[2] += 1

    variant.sequence = choose_best_seq(reference_sequences)
    if not variant.sequence:
        variant.correct = False
        print('Skipping: No sequence for ', variant.refsnp_id, 'variant.')

    ref = variant.sequence[offset:-offset + variant.length]
    # print('Context: ' + show_pos_with_context(seq, offset, -offset))

    vcf_data = get_vcf_by_variant(pos, variant)

    if not vcf_data:
        print('Skipping: Lack of VCF data for', variant.refsnp_id, 'variant.')
        variant.correct = False
        return False

    assert len(vcf_data.ALT) == 1

    variant.alt = str(vcf_data.ALT[0])
    variant.ref = str(vcf_data.REF)
    variant.gene = vcf_data.INFO['GENE'][0]

    if variant.ref != ref:
        print(
            'VCF says ref is', variant.ref, 'but analysis pointed to',
            ref, 'for', variant.refsnp_id
        )

    if variant.refsnp_source != 'COSMIC':
        # Allele is not informative for entries from cosmic ('COSMIC_MUTATION')
        if variant.ref != variant.allele_1:
            print(
                'VCF says ref is', variant.ref, 'but biomart believes it\'s',
                variant.allele_q, 'for', variant.refsnp_id
            )

    analyze_poly_a(variant, offset)

    return True


def analyze_poly_a(variant, offset):

    ref_seq = variant.sequence
    alt = variant.alt
    mutated_seq = ref_seq[:offset] + str(alt) + ref_seq[offset + variant.length:]
    o.print('Mutated: ' + show_pos_with_context(mutated_seq, offset, -offset))

    has_aaa, before_len = poly_a(
        ref_seq,
        offset,
        len(ref_seq) - offset
    )

    will_have, after_len = poly_a(
        mutated_seq,
        offset,
        len(mutated_seq) - offset
    )

    variant.has_poly_a = has_aaa
    variant.will_have_poly_a = will_have
    variant.poly_aaa_before = before_len
    variant.poly_aaa_after = after_len
    variant.poly_aaa_change = variant.poly_aaa_after - variant.poly_aaa_before
    variant.poly_aaa_increase = variant.poly_aaa_after > variant.poly_aaa_before
    variant.poly_aaa_decrease = variant.poly_aaa_after < variant.poly_aaa_before


genes_from_patacsdb = gene_names_from_patacsdb_csv(None)

dataset = BiomartDataset(
    'http://www.ensembl.org/biomart',
    name='hsapiens_snp_som'
)


@cached(action='load')
def cachable_variants_by_gene():
    return get_variants_by_genes(dataset, genes_from_patacsdb)


variants_by_gene = cachable_variants_by_gene()


@cached(action='load')
def cachable_transcripts_to_load():
    return get_all_used_transcript_ids(variants_by_gene)


transcripts_to_load = cachable_transcripts_to_load()


@cached(action='load')
def cachable_cds_db():
    return SequenceDB(
        index_by='transcript',
        sequence_type='cds',
        restrict_to=transcripts_to_load
    )


@cached(action='load')
def cachable_cdna_db():
    return SequenceDB(
        index_by='transcript',
        sequence_type='cdna',
        restrict_to=transcripts_to_load
    )


@cached(action='load')
def cachable_dna_db():
    chromosomes = map(str, range(1, 23)) + ['X', 'Y', 'MT']
    dna_db = {}

    for chromosome in chromosomes:
        dna_db[chromosome] = FastSequenceDB(
            sequence_type='dna',
            id_type='chromosome.' + chromosome
        )
    return dna_db


def analyze_variant_here(variant):

    dna_db = cachable_dna_db()
    cds_db = cachable_cds_db()
    cdna_db = cachable_cdna_db()

    try:
        analyze_variant(
            variant,
            cds_db,
            cdna_db,
            dna_db
            # vcf_cosmic,
            # vcf_ensembl,
            # hgnc_by_ensembl
        )

    except Exception:
        traceback.print_exc()
        variant.correct = False

    del dna_db, cds_db, cdna_db

    return variant


def parse_variants(variants_by_gene):
    from multiprocessing import Pool
    """
    vcf.parser uses 0-based coordinates:
    http://pyvcf.readthedocs.org/en/latest/_modules/vcf/parser.html?highlight=coordinates
    ensembl uses 1-based coordinates:
    http://www.ensembl.org/info/docs/api/core/core_tutorial.html#coordinates
    """

    variants_by_gene_by_transcript = {}

    all_variants_count = 0

    cosmic_genes_to_load = set()

    parsing_pool = Pool(maxtasksperchild=1)
    print('Parsing variants:')
    for gene, variants in tqdm(variants_by_gene.iteritems(), total=len(variants_by_gene)):

        # Just to be certain
        variants_unique_ids = set(variant.refsnp_id for variant in variants)
        if len(variants_unique_ids) != len(variants):
            raise Exception('Less unique ids than variants!')

        all_variants_count += len(variants)

        variants = parsing_pool.map(analyze_variant_here, variants)
        """
        for variant in variants:
            analyze_variant(
                variant,
                cds_db,
                cdna_db,
                dna_db,
                vcf_cosmic,
                vcf_ensembl
            )
        """
        gc.collect()

        # Remove variants with non-complete data
        correct_variants = filter(lambda variant: variant.correct, variants)

        # Remove variants from source other than COSMIC (later I am processing
        # cosmic-specific data so db_snp variants are not relevant here)
        correct_variants = filter(
            lambda variant: variant.refsnp_source == 'COSMIC',
            list(correct_variants)
        )

        cosmic_genes_to_load.update(
            [
                variant.vcf_gene
                for variant in correct_variants
            ]
        )

        by_transcript = {}
        # czwartki koło południa CENT

        for variant in correct_variants:
            # The problem with the ensembl's biomart is that it returns records
            # from cosmic without information about the transcript, so we have
            # often a few identical records with only the refsnp_id different,
            # as for example: COSM3391893, COSM3391894
            # Fortunately the transcript id is encoded inside vcf_data retrived
            # from biomart inside the gene identifer (if it is abset, then we
            # have a canonical transcript, at least it is the best guess), eg.:
            # ANKRD26_ENST00000376070 | ANKRD26_ENST00000436985 | ANKRD26

            gene_transcript_id = variant.vcf_gene
            transcript_id = gene_transcript_id

            try:
                by_transcript[transcript_id].append(variant)
            except KeyError:
                by_transcript[transcript_id] = [variant]

        variants_by_gene_by_transcript[gene] = by_transcript

    o.print('All variants', all_variants_count)

    parsing_pool.close()

    return variants_by_gene_by_transcript, cosmic_genes_to_load


def report(name, data, comment=None):
    """Generates list-based report quickly.

    File will be placed in 'reports' dir and name will be derived from 'name' of
    the report. The data should be an iterable collection of strings.
    Empty reports will not be created.
    """
    directory = 'reports_three'

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


def get_all_variants(variants_by_transcript, gene):

    unique_variants = {}
    duplicated = set()

    for gene_transcript_id, variants in variants_by_transcript.iteritems():
        for variant in variants:
            key = '\t'.join([
                variant.chr_name,
                variant.chrom_start,
                variant.chrom_end,
                # variant.ref, # TODO
                str(variant.vcf_data.REF),
                ','.join([str(n) for n in variant.alt]),  # just variant TODO
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

    report('Duplicated records for gene: ' + gene, duplicated)

    return unique_variants.values()


def select_poly_a_related_variants(variants):
    return [
        variant
        for variant in variants
        if variant.has_poly_a or variant.will_have_poly_a
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
                variant.poly_aaa_increase,
                variant.poly_aaa_decrease,
                variant.poly_aaa_change
            ]))
            for variant in poly_a_related_variants
        ]

    report(
        'poly aaa increase and decrease by variants',
        variant_aaa_report,
        'snp_id\tpoly_aaa_increase\tpoly_aaa_decrease\tpoly_aaa_change'
    )


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
        variants = get_all_variants(variants_by_transcript, gene)

        variants_count += len(variants)

        o.print(gene, 'with its', len(variants), 'variants analysed')

        poly_a_variants = filter(
            lambda variant: variant.has_poly_a,
            variants
        )

        poly_a_potential_variants = filter(
            lambda variant: variant.will_have_poly_a,
            variants
        )

        # loss of poly_aaa, decrease in length (-1 per residue)
        decrease = 0
        # gain of poly_aaa: , increase in length (+1 per residue)
        increase = 0

        # the same as: poly_a_related_variants = select_poly_a_related_variants(variants)
        poly_a_related_variants = list(poly_a_variants + poly_a_potential_variants)
        poly_a_related_variants_count += len(poly_a_related_variants)

        # give scores for length increase
        for variant in poly_a_related_variants:
            if variant.poly_aaa_increase:
                increase += 1
            elif variant.poly_aaa_decrease:
                decrease += 1

        assert increase >= len(poly_a_potential_variants) - len(poly_a_variants)

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


def poly_aaa_vs_expression(variants_by_gene_by_transcript, cache_action='load'):

    from expression_database import ExpressionDatabase
    from expression_database import import_expression_data

    bdb = ExpressionDatabase('expression_slope_in_tissues_by_mutation.db')

    if cache_action == 'save':
        import_expression_data(bdb)

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
        variants = get_all_variants(variants_by_transcript, gene)

        poly_a_related_variants = select_poly_a_related_variants(variants)

        for variant in poly_a_related_variants:

            expression_data = bdb.get_by_mutation(variant)

            if not expression_data:
                # print('No expression for', variant.refsnp_id)
                variant.expression_trend = 'none'
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
                    variant.expression_trend = 'up'
                else:
                    variant.expression_trend = 'rather_up'
            # is this rather down?
            elif len(expression_down) > len(expression_up):
                # is this certainly down?
                if is_length_difference_big(expression_down, expression_up):
                    variant.expression_trend = 'down'
                else:
                    variant.expression_trend = 'rather_down'
            # is unchanged?
            else:
                variant.expression_trend = 'constant'

            variant.expression_up_in_X_cases = len(expression_up)
            variant.expression_down_in_X_cases = len(expression_down)

            gtex_report += [(
                variant.refsnp_id,
                variant.expression_up_in_X_cases,
                variant.expression_down_in_X_cases,
                variant.expression_trend,
                variant.poly_aaa_increase,
                variant.poly_aaa_decrease,
                variant.poly_aaa_change
            )]

        gtex_report_by_genes += [(
            gene,
            sum('up' in v.expression_trend for v in poly_a_related_variants),
            sum('down' in v.expression_trend for v in poly_a_related_variants),
            sum(v.poly_aaa_increase for v in poly_a_related_variants),
            sum(v.poly_aaa_decrease for v in poly_a_related_variants)
        )]

    report(
        'Expression table for variants (based on data from GTEx)',
        ['\t'.join(map(str, line)) for line in gtex_report],
        'variant\texpression+\texpression-\ttrend\taaa+\taaa-\taaa change'
    )

    report(
        'Expression table for genes (based on data from GTEx)',
        ['\t'.join(map(str, line)) for line in gtex_report_by_genes],
        'gene\texpression+\texpression-\t#aaa+\t#aaa-'
    )

    print('Done')


def main(args, dataset):
    """The main workflow happens here."""
    cache = args.cache

    genes_from_patacsdb = gene_names_from_patacsdb_csv(args.number)

    @cached(action='load')
    def cachable_variants_by_gene():
        return get_variants_by_genes(dataset, genes_from_patacsdb)

    variants_by_gene = cachable_variants_by_gene()

    @cached(action='load')
    def cachable_transcripts_to_load():
        return get_all_used_transcript_ids(variants_by_gene)

    transcripts_to_load = cachable_transcripts_to_load()

    @cached(action='load')
    def cachable_cds_db():
        return SequenceDB(
            index_by='transcript',
            sequence_type='cds',
            restrict_to=transcripts_to_load
        )

    #cds_db = cachable_cds_db()

    @cached(action='load')
    def cachable_cdna_db():
        return SequenceDB(
            index_by='transcript',
            sequence_type='cdna',
            restrict_to=transcripts_to_load
        )

    import traceback
    import sys
    import cPickle as pickle

    #cdna_db = cachable_cdna_db()

    if cache == 'save':

        variants_by_gene_by_transcript, cosmic_genes_to_load = parse_variants(
            variants_by_gene
        )

        try:
            with open('.variants_by_gene_by_transcript.cache', 'wb') as f:
                pickle.dump(variants_by_gene_by_transcript, f, protocol=pickle.HIGHEST_PROTOCOL)
            with open('.cosmic_genes_to_load.cache', 'wb') as f:
                pickle.dump(cosmic_genes_to_load, f, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception:
            traceback.print_exc()

    else:
        with open('.variants_by_gene_by_transcript.cache', 'rb') as f:
            variants_by_gene_by_transcript = pickle.load(f)
        with open('.cosmic_genes_to_load.cache', 'rb') as f:
            cosmic_genes_to_load = pickle.load(f)

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
                    'cosmic/CosmicCompleteCNA.tsv',
                    restrict_to=cosmic_genes_to_load
                )

            cna = cachable_cna()
            summarize_copy_number_expression(variants_by_gene_by_transcript, cna)
    except Exception:
        traceback.print_exc()

    try:
        if 'poly_aaa_vs_expression' in args.report:
            poly_aaa_vs_expression(variants_by_gene_by_transcript, cache)
    except Exception:
        traceback.print_exc()


if __name__ == '__main__':

    import argparse

    to_show = ['databases', 'datasets', 'filters', 'attributes',
               'attributes_by_page', 'some_variant']
    cache_actions = ['load', 'save']
    parser = argparse.ArgumentParser(description='Find SNPs')
    parser.add_argument('--show', choices=to_show)
    parser.add_argument('--cache', choices=cache_actions)
    parser.add_argument('--report', nargs='+', choices=[
        'list_poly_aaa_variants', 'copy_number_expression',
        'poly_aaa_vs_expression'
    ])
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

    subcommands = ['show', 'cache', 'profile']
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
            main(args, snp_dataset)
