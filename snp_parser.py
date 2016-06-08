#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import vcf
import sys
from poly_a import poly_a
from fasta_sequence_db import SequenceDB, FastSequenceDB
from output_formatter import OutputFormatter
from biomart_data import BiomartData, BiomartDataset
from cna_by_transcript import CompleteCNA

# pip install pyvcf pysam biomart

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
            # 'validated',
            # 'minor_allele_freq',
            # 'minor_allele_count',
            # 'variation_names',
            'phenotype_description',
            'set_name',
            'set_description',
            # 'clinical_significance',
            # 'refsnp_source_description',
            'cdna_start',
            'cdna_end',
            # 'consequence_type_tv',
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
                    'synonymous_variant',
                    'missense_variant',
                    'stop_gained',
                    'coding_sequence_variant'
                ]
        })
        super(self.__class__, self).__init__(dataset, attributes, filters)


def show_pos_with_context(seq, start, end):
    return seq[:start] + '→' + seq[start:end] + '←' + seq[end:]


def gene_names_from_patacsdb_csv(how_many):

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


def parse(dataset, how_many_genes=False):

    genes_from_patacsdb = gene_names_from_patacsdb_csv(how_many_genes)

    variants_by_gene = {}

    previous_id = None

    for start in range(0, len(genes_from_patacsdb), 300):

        filters = {'ensembl_gene': genes_from_patacsdb[start:start + 300]}

        variants = VariantsData(filters=filters, dataset=dataset)

        # Those attributes are only part that changes acros 'duplicated'
        # records returned by ensembl's biomart for a given mutation.
        # To avoid redundancy I am assembling all redundant entries into
        # single one, grouped by refsnp_id (in case of cosmic mutations,
        # the cosmic id is stored inside this field). Attributes that are
        # chaning acros 'duplicated' entries are sotred in a list so no
        # information is lost. Nontheless this does not reselve all the
        # issues with redundancy of retrived data - some are processed later.
        redundant_attrs = ['phenotype_description', 'set_name', 'set_description']

        # I assert that two redundant entires comes always one after another.
        for variant in variants:

            gene = variant.ensembl_gene_stable_id

            assert variant.refsnp_id

            if variant.refsnp_id == previous_id:
                last = variants_by_gene[gene][-1]
                assert previous_id == last.refsnp_id
                for attr in redundant_attrs:
                    last.get(attr).append(variant.get(attr))
            else:
                for attr in redundant_attrs:
                    variant.set(attr, [variant.get(attr)])
                if gene not in variants_by_gene:
                    variants_by_gene[gene] = []
                variants_by_gene[gene].append(variant)
                previous_id = variant.refsnp_id

        variants.iterator = None

    return variants_by_gene


def ref_seq_len(src, ref):
    if src not in ref:
        return 0
    return len(ref[src].strip('-'))


def inspect(obj):
    """
    Just for debugging and exploration
    """
    o.unmute()
    o.print(obj)
    o.indent()
    o.print(type(obj))
    for k, v in obj.__dict__.items():
        o.print(k, v)
    o.print(dir(obj))
    o.outdent()


def get_vcf_by_variant(vcf, pos, variant):

    if variant.refsnp_source == 'COSMIC':
        record_id = variant.refsnp_id
    else:

        data = BiomartData(
            dataset=BiomartDataset(BIOMART_URL, name='hsapiens_gene_ensembl'),
            attributes=['hgnc_symbol'],
            filters={
                'ensembl_transcript_id': variant.ensembl_transcript_stable_id
            }
        )

        name = list(data)[0].hgnc_symbol

        allele = variant.allele.replace('/', '>')
        record_id = ''.join([name, ':c.', variant.cds_start, allele])

    return get_vcf_by_id(vcf, pos, record_id)


def get_vcf_by_id(vcf, pos, record_id):

    vcf_data = None

    for record in vcf.fetch(*pos):
        if record.ID == record_id:
            vcf_data = record
            break

    return vcf_data


def analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl):

    offset = 20

    o.mute()
    o.print('Variant name:', variant.refsnp_id)

    transcript_id = variant.ensembl_transcript_stable_id
    strand = int(variant.ensembl_transcript_chrom_strand)

    if variant.refsnp_source != 'COSMIC':
        print('Found variant', variant.refsnp_id,
              'from source other than COSMIC:', variant.refsnp_source)

    reference_nuc = {}
    reference_seq = {}

    for db in [cdna_db, cds_db]:
        src = db.sequence_type

        try:
            start, end = variant.get(src + '_start'), variant.get(src + '_end')
        except IndexError:
            o.print('Lack of', src, 'coordinates for variant:',
                    variant.refsnp_id, 'in context of',
                    transcript_id)
            continue

        seq = db.fetch(transcript_id, strand, start, end, offset)

        if not db.has(transcript_id) or not seq:
            o.print('Lack of transcript in', src)
            continue

        reference_nuc[src] = seq[offset:-offset]
        reference_seq[src] = seq

    # o.print('Ancestral allele: ' + variant.allele_1)
    # Allele is not informative for entries from cosmic ('COSMIC_MUTATION')
    # o.print('Allele: ' + variant.allele)
    # o.print('consequence_allele_string:' + variant.consequence_allele_string)

    reference_nuc['biomart (ancestral)'] = variant.allele_1

    chromosome = dna_db[variant.chr_name]

    start, end = chromosome.parse_coordinates(variant.chrom_start, variant.chrom_end)

    pos = [str(variant.chr_name), int(variant.chrom_start), int(variant.chrom_end)]

    seq = chromosome.fetch(pos[1], pos[2], offset)
    # nuc = chromosome.fetch(pos[1], pos[2])

    reference_nuc['genome'] = seq[offset:-offset]
    reference_seq['genome'] = seq

    # o.print('Chromosomal position: {0}:{1}-{2}'.format(*pos))

    # to get to vcf stored data by vcf reader, change coordinates to 0-based
    pos[1] -= 1
    pos[2] -= 1

    # and represent them as range (if you had n:n pointing to a single base, use n:n+1)
    pos[2] += 1

    for src, seq in reference_seq.items():
        if not seq:
            del reference_seq[src]
            # o.print('Lack of reference sequence in ' + src)

    temp_ref_seq = reference_seq['genome']
    consistent = True

    for src, seq in reference_seq.items():
        if seq.startswith('-') or seq.endswith('-'):
            o.print('Offset surpasses {0} transcript span'.format(src))
        while seq.startswith('-'):
            seq = seq[1:]
            temp_ref_seq = temp_ref_seq[1:]
        while seq.endswith('-'):
            seq = seq[:-1]
            temp_ref_seq = temp_ref_seq[:-1]
        if seq != temp_ref_seq:
            consistent = False

    if not consistent:
        o.print('Reference sequences are not consistent:')
        for src, seq in reference_seq.items():
            if seq:
                o.print(src, ':\t', show_pos_with_context(seq, offset, -offset))

    variant.cds_cdna_inconsistent = False
    if reference_seq.get('cds', '') != reference_seq.get('cdna', ''):
        cdna_real_len = ref_seq_len('cdna', reference_seq)
        cds_real_len = ref_seq_len('cds', reference_seq)
        consensus = cds_real_len != cdna_real_len
        if 'cds' in reference_seq and 'cdna' in reference_seq:
            cds = reference_seq['cds']
            cdna = reference_seq['cds']
            consensus = True
            if len(cdna) == len(cds):
                for i in range(len(cdna)):
                    if cdna[i] != '-' and cds[i] != '-' and cds[i] != cdna[i]:
                        consensus = False
                        break
        if not consensus:
            o.unmute()
            o.print(reference_seq)
            o.print('cdna and cds sequences are totally inconsistent')
            exit()
        else:
            o.print('cds and cdna of different length')
            variant.cds_cdna_inconsistent = True

    if (ref_seq_len('cds', reference_seq) >= ref_seq_len('cdna', reference_seq) and
            ref_seq_len('cds', reference_seq)):
        ref_seq = reference_seq['cds']
        o.print('Chosing cds sequence as reference')
    elif ref_seq_len('cdna', reference_seq):
        ref_seq = reference_seq['cdna']
        o.print('Chosing cdna sequence as reference')
    else:
        ref_seq = reference_seq['genome']
        o.print('Chosing genome sequence as reference')

    variant.sequence = ref_seq
    o.print('Context: ' + show_pos_with_context(ref_seq, offset, -offset))

    main_vcf = vcf_cosmic if variant.refsnp_source == 'COSMIC' else vcf_ensembl

    vcf_data = get_vcf_by_variant(main_vcf, pos, variant)

    if not vcf_data:
        print('Lack of VCF data for', variant.refsnp_id, 'variant. Skipping')
        variant.correct = False
        return False
    else:
        variant.correct = True

    assert len(vcf_data.ALT) == 1
    variant.vcf_data = vcf_data

    variant.alt = str(vcf_data.ALT)

    analyze_poly_a(variant, offset)

    o.unmute()

    return True


def analyze_poly_a(variant, offset):

    ref_seq = variant.sequence
    alt = variant.alt
    mutated_seq = ref_seq[:offset] + str(alt) + ref_seq[offset + 1:]
    o.print('Mutated: ' + show_pos_with_context(mutated_seq, offset, -offset))

    has_aaa, before_len = poly_a(ref_seq, offset, len(ref_seq) - offset)
    will_have, after_len = poly_a(mutated_seq, offset, len(mutated_seq) - offset)

    variant.has_poly_a = has_aaa
    variant.will_have_poly_a = will_have
    variant.poly_aaa_before = before_len
    variant.poly_aaa_after = after_len


def parse_variants(cds_db, cdna_db, variants_by_gene):

    """
    vcf.parser uses 0-based coordinates:
    http://pyvcf.readthedocs.org/en/latest/_modules/vcf/parser.html?highlight=coordinates
    ensembl uses 1-based coordinates:
    http://www.ensembl.org/info/docs/api/core/core_tutorial.html#coordinates
    """

    chromosomes = map(str, range(1, 23)) + ['X', 'Y', 'MT']

    dna_db = {}
    for chromosome in chromosomes:
        dna_db[chromosome] = FastSequenceDB(sequence_type='dna', id_type='chromosome.' + chromosome)

    """
    Przygotowanie pliku tabix:
    1. zainstaluj htslib z http://www.htslib.org
    2. run ./create_tabix.sh filename
    """
    vcf_ensembl = vcf.Reader(filename='ensembl_vcf/Homo_sapiens.vcf.gz')
    vcf_cosmic = vcf.Reader(filename='cosmic/CosmicCodingMuts.vcf.gz')

    variants_by_gene_by_transcript = {}

    all_variants_count = 0

    cosmic_genes_to_load = set()

    for gene, variants in variants_by_gene.iteritems():

        # Just to be certain
        variants_unique_ids = set(variant.refsnp_id for variant in variants)
        if len(variants_unique_ids) != len(variants):
            raise Exception('Less unique ids than variants!')

        all_variants_count += len(variants)

        for variant in variants:
            analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl)

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
                variant.vcf_data.INFO['GENE']
                for variant in correct_variants
            ]
        )

        by_transcript = {}

        for variant in correct_variants:
            # The problem with the ensembl's biomart is that it returns records
            # from cosmic without information about the transcript, so we have
            # often a few identical records with only the refsnp_id different,
            # as for example: COSM3391893, COSM3391894
            # Fortunately the transcript id is encoded inside vcf_data retrived
            # from biomart inside the gene identifer (if it is abset, then we
            # have a canonical transcript, at least it is the best guess), eg.:
            # ANKRD26_ENST00000376070 | ANKRD26_ENST00000436985 | ANKRD26

            gene_transcript_id = variant.vcf_data.INFO['GENE']
            transcript_id = gene_transcript_id

            try:
                by_transcript[transcript_id].append(variant)
            except KeyError:
                by_transcript[transcript_id] = [variant]

        variants_by_gene_by_transcript[gene] = by_transcript

    o.print('All variants', all_variants_count)

    return variants_by_gene_by_transcript, cosmic_genes_to_load


def report(name, data, comment=None):
    if not data:
        return
    with open('reports/' + name.replace(' ', '_') + '.txt', 'w') as f:
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
                variant.vcf_data.REF,
                ','.join([str(n) for n in variant.vcf_data.ALT]),
                variant.ensembl_gene_stable_id,
                gene_transcript_id
            ])
            if key not in unique_variants.keys():
                variant.affected_transcripts = set(gene_transcript_id)
                unique_variants[key] = variant
            else:
                unique_variants[key].affected_transcripts.add(gene_transcript_id)
                stored_variant = unique_variants[key]
                duplicated.update([variant.refsnp_id], [stored_variant.refsnp_id])

    report('Duplicated records for gene: ' + gene, duplicated)

    return unique_variants.values()


def summarize(variants_by_gene_by_transcript, cna):

    variants_count = 0
    poly_a_related_variants_count = 0

    to_report = []
    import operator

    no_expression = set()

    for gene, variants_by_transcript in variants_by_gene_by_transcript.iteritems():

        expression = [0, 0, 0]
        for gene_transcript_id, variants in variants_by_transcript.iteritems():

            try:
                expr = cna.get_by_gene_and_transcript(gene_transcript_id)
                expression = map(operator.add, expression, expr)
            except KeyError:
                # no expression data for this transcript
                no_expression.add(gene_transcript_id)
                continue

        # Treat all variants the same way - just remove duplicates
        variants = get_all_variants(variants_by_transcript, gene)

        variants_count += len(variants)

        o.print(gene, 'with its', len(variants), 'variants analysed')

        poly_a_variants = filter(lambda variant: variant.has_poly_a, variants)
        poly_a_potential_variants = filter(lambda variant: variant.will_have_poly_a, variants)

        # loss of poly_aaa, decrease in length (-1 per residue)
        decrease = 0
        # gain of poly_aaa: , increase in length (+1 per residue)
        increase = 0

        poly_a_related_variants = list(poly_a_variants + poly_a_potential_variants)
        poly_a_related_variants_count += len(poly_a_related_variants)

        # give scores for length increase
        for variant in poly_a_related_variants:
            if variant.poly_aaa_after > variant.poly_aaa_before:
                increase += 1
            elif variant.poly_aaa_after < variant.poly_aaa_before:
                decrease += 1

        assert increase >= len(poly_a_potential_variants) - len(poly_a_variants)

        to_report += [(gene, expression[0], expression[2], increase, decrease)]

    report('No expression data for some transcripts', no_expression)

    gene_count = len(variants_by_gene_by_transcript)

    print('Analyzed genes:', gene_count)
    print('Analyzed variants:', variants_count)
    print('poly_a related variants', poly_a_related_variants_count)

    report('Poly A and expression table',
           ['\t'.join(map(str, line)) for line in to_report],
           'Gene\tCNV+\tCNV-\tAAA+\tAAA-')


def get_unique_transcript_ids(variants_by_gene):

    transcripts_to_load = []
    for gene, variants in variants_by_gene.iteritems():
        transcripts_to_load += [v.ensembl_transcript_stable_id for v in variants]

    return set(transcripts_to_load)


def parse_dataset(dataset):
    o.print('Parsing dataset:', dataset.name)
    variants_by_gene = parse(dataset)

    transcripts_to_load = get_unique_transcript_ids(variants_by_gene)

    cds_db = SequenceDB(
        index_by='transcript',
        sequence_type='cds',
        restrict_to=transcripts_to_load)
    cdna_db = SequenceDB(
        index_by='transcript',
        sequence_type='cdna',
        restrict_to=transcripts_to_load)

    o.print('Dataset parsed')

    return variants_by_gene, cds_db, cdna_db


def main(args, dataset):
    import cPickle as pickle

    cache = args.cache
    cache_name = '.cache'

    if cache == 'load':
        with open(cache_name, 'rb') as f:
            parsed_data = pickle.load(f)
        o.print('Variants data loaded from cache')
    else:
        parsed_data = parse_dataset(dataset)

        if cache == 'save':
            with open(cache_name, 'wb') as f:
                pickle.dump(parsed_data, f, protocol=pickle.HIGHEST_PROTOCOL)
            o.print('variants data saved to cache')

    variants_by_gene, cds_db, cdna_db = parsed_data

    variants_by_gene_by_transcript, cosmic_genes_to_load = parse_variants(cds_db, cdna_db, variants_by_gene)

    if cache == 'load':
        with open(cache_name + '-cna', 'rb') as f:
            cna = pickle.load(f)
        o.print('CNA loaded from cache')
    else:
        cna = CompleteCNA(
            'cosmic/CosmicCompleteCNA.tsv',
            restrict_to=cosmic_genes_to_load)
        if cache == 'save':
            with open(cache_name + '-cna', 'wb') as f:
                pickle.dump(cna, f, protocol=pickle.HIGHEST_PROTOCOL)

    summarize(variants_by_gene_by_transcript, cna)


if __name__ == '__main__':

    import argparse

    to_show = ['databases', 'datasets', 'filters', 'attributes',
               'attributes_by_page', 'some_variant']
    cache_actions = ['load', 'save']
    parser = argparse.ArgumentParser(description='Find SNPs')
    parser.add_argument('--show', choices=to_show)
    parser.add_argument('--cache', choices=cache_actions)
    parser.add_argument('--profile', action='store_true')
    parser.add_argument(
        '--dataset',
        type=str,
        help='name of biomart dataset to be used, eg. hsapiens_snp',
        default='hsapiens_snp_som')
    parser.add_argument(
        '--biomart',
        type=str,
        help='URL of biomart to be used. '
             'For ensembl mirrors replace www with: uswest, useast or asia',
        default='http://www.ensembl.org/biomart')

    subcommands = ['show', 'cache', 'profile']
    arguments = ['--' + a if a in subcommands else a for a in sys.argv[1:]]

    args = parser.parse_args(arguments)

    global BIOMART_URL
    BIOMART_URL = args.biomart

    snp_dataset = BiomartDataset(args.biomart, name=args.dataset)

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
            variants_by_gene = parse(snp_dataset, how_many_genes=5)
            gene, variants = variants_by_gene.popitem()
            print(variants[0])
    else:

        if args.profile:
            import profile
            profile.run('main(args, snp_dataset)')
        else:
            main(args, snp_dataset)
