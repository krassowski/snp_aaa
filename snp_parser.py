#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

DATASET = 'hsapiens_snp_som'
BIOMART_URL = 'http://www.ensembl.org/biomart'
"""
Mirrors:
    uswest.ensembl.org
    useast.ensembl.org
    asia.ensembl.org
"""

# pip install pyvcf pysam biomart
import vcf
import sys

from poly_a import has_poly_a
from fasta_sequence_db import SequenceDB, FastSequenceDB
from output_formatter import OutputFormatter
from copy_number_map import CopyNumberMap
from biomart_data import BiomartData, BiomartDataset

o = OutputFormatter()


def inspect(obj):
    """
    Just for debugging and exploration
    """
    o.print(obj)
    o.indent()
    o.print(type(obj))
    for k, v in obj.__dict__.items():
        o.print(k, v)
    o.print(dir(obj))
    o.outdent()


class VariantsData(BiomartData):

    def __init__(self, dataset=None, attributes=None, filters=None):

        if not dataset:
            dataset = BiomartDataset(BIOMART_URL, name=DATASET)

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
            'validated',
            'minor_allele_freq',
            'minor_allele_count',
            'variation_names',
            'phenotype_description',
            'set_name',
            'set_description',
            'clinical_significance',
            'refsnp_source_description',
            'cdna_start',
            'cdna_end',
            'consequence_type_tv',
            'ensembl_gene_stable_id',
            'ensembl_transcript_stable_id',
            'ensembl_transcript_chrom_strand',
            'cds_start',
            'cds_end',
            'consequence_type_tv',
            'consequence_allele_string'
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
    genes = []
    count = 0
    with open('patacsdb_all.csv', 'r') as f:
        for line in f:
            if count >= how_many:
                break
            count += 1
            if line:
                data = [x.strip() for x in line.split(',')]
                genes.append(data[1])
    return genes


def parse(how_many=250):

    genes_from_patacsdb = gene_names_from_patacsdb_csv(how_many)

    filters = {
        'ensembl_gene': genes_from_patacsdb
    }

    variants = VariantsData(filters=filters)

    variants_by_gene = {}

    previous_id = None

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

    return variants_by_gene, variants


def ref_seq_len(src, ref):
    if src not in ref:
        return 0
    return len(ref[src].strip('-'))


def analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl):
    o.mute()

    o.print('Variant name:', variant.refsnp_id)
    o.indent()

    # gene_id = variant.ensembl_gene_stable_id

    # trasnscript is a transcript of a gene (reference), not of the variant!
    transcript_id = variant.ensembl_transcript_stable_id
    strand = int(variant.ensembl_transcript_chrom_strand)

    offset = 20

    o.print('Transcript', transcript_id)
    o.indent()

    if variant.refsnp_source != 'COSMIC':
        o.unmute()
    o.print(variant.consequence_type_tv)
    o.print('strand:', strand)

    reference_nuc = {}
    reference_seq = {}

    for db in [cdna_db, cds_db]:
        src = db.sequence_type

        start, end = variant.get(src + '_start'), variant.get(src + '_end')
        seq = db.fetch(transcript_id, strand, start, end, offset)

        if not db.has(transcript_id) or not seq:
            o.print('Lack of transcript in', src)
            continue

        reference_nuc[src] = seq[offset:-offset]
        reference_seq[src] = seq

    o.outdent()
    o.print('Ancestral allele: ' + variant.allele_1)
    # Allele is usually not informative for cosmic sourced entries ('COSMIC_MUTATION')
    o.print('Allele: ' + variant.allele)
    o.print('consequence_allele_string:' + variant.consequence_allele_string)

    reference_nuc['biomart (ancestral)'] = variant.allele_1

    chromosome = dna_db[variant.chr_name]

    start, end = chromosome.parse_coordinates(variant.chrom_start, variant.chrom_end)

    pos = [str(variant.chr_name), int(variant.chrom_start), int(variant.chrom_end)]

    seq = chromosome.fetch(pos[1], pos[2], offset)
    # nuc = chromosome.fetch(pos[1], pos[2])

    reference_nuc['genome'] = seq[offset:-offset]
    reference_seq['genome'] = seq

    o.print('Chromosomal position: {0}:{1}-{2}'.format(*pos))

    # to get to vcf stored data by vcf reader, change coordinates to 0-based
    pos[1] -= 1
    pos[2] -= 1

    # and represent them as range (eg. if you had n:n pointing to a single base, use n:n+1)
    pos[2] += 1

    for src, seq in reference_seq.items():
        if not seq:
            del reference_seq[src]
            o.print('Lack of reference sequence in ' + src)

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
        # TODO: skip?
        o.print('Reference sequences are not consistent:')
        for src, seq in reference_seq.items():
            if seq:
                o.print(src, ':\t', show_pos_with_context(seq, offset, -offset))

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

    if ref_seq_len('cds', reference_seq) >= ref_seq_len('cdna', reference_seq) and ref_seq_len('cds', reference_seq):
        ref_seq = reference_seq['cds']
        o.print('Chosing cds sequence as reference')
    elif ref_seq_len('cdna', reference_seq):
        ref_seq = reference_seq['cdna']
        o.print('Chosing cdna sequence as reference')
    else:
        ref_seq = reference_seq['genome']
        o.print('Chosing genome sequence as reference')

    o.print('Context: ' + show_pos_with_context(ref_seq, offset, -offset))

    main_vcf = vcf_cosmic if variant.refsnp_source == 'COSMIC' else vcf_ensembl

    vcf_data = None

    for record in main_vcf.fetch(*pos):

        if record.ID == variant.refsnp_id:
            vcf_data = record
            break

    if not vcf_data:
        o.unmute()
        o.print('Lack of VCF data for', variant.refsnp_id, 'variant. Skipping')
        # eg: cancer.sanger.ac.uk/cosmic/mutation/overview?id=4967608
        variant.correct = False

        o.outdent()
        return False

    # o.print(record.num_called, record.call_rate, record.num_unknown)  # ZeroDivisionError ?
    # o.print(record.num_hom_ref, record.num_het, record.num_hom_alt)
    # o.print(record.nucl_diversity, record.aaf, record.heterozygosity)  # ZeroDivisionError ?
    # o.print(record.is_snp, record.is_indel, record.is_transition, record.is_deletion)
    # o.print(record.is_monomorphic)
    # o.print(vcf_data)

    alt = vcf_data.ALT

    assert len(alt) == 1
    alt = alt[0]

    mutated_seq = ref_seq[:offset] + str(alt) + ref_seq[offset + 1:]

    variant.sequence = ref_seq
    variant.vcf_data = vcf_data
    variant.has_poly_a = has_poly_a(ref_seq, offset, len(ref_seq) - offset)
    variant.will_have_poly_a = has_poly_a(mutated_seq, offset, len(mutated_seq) - offset)
    variant.poly_a_decrease = variant.has_poly_a > variant.will_have_poly_a
    variant.poly_a_increase = variant.has_poly_a < variant.will_have_poly_a
    variant.correct = True

    o.print('Mutated: ' + show_pos_with_context(mutated_seq, offset, -offset))

    o.unmute()
    o.outdent()

    if variant.allele != 'COSMIC_MUTATION':
        exit()

    return True


def parse_variants(cds_db, cdna_db):

    """
    vcf.parser uses 0-based coordinates:
    http://pyvcf.readthedocs.org/en/latest/_modules/vcf/parser.html?highlight=coordinates
    ensembl uses 1-based coordinates:
    http://www.ensembl.org/info/docs/api/core/core_tutorial.html#coordinates
    """

    chromosomes = map(str, range(1, 23)) + ['X', 'Y']

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
        assert len(variants_unique_ids) == len(variants)

        all_variants_count += len(variants)

        for variant in variants:
            analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl)

        # Remove variants with non-complete data 
        correct_variants = filter(lambda variant: variant.correct, variants)

        cosmic_genes_to_load.update([variant.vcf_data.INFO['GENE'] for variant in correct_variants])

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
            """
            chunks = gene_transcript_id.split('_')
            if len(chunks) > 1:
                transcript_id = chunks[1]
                # just in case to detect potential ambuguity
                assert len(chunks) == 2
            else:
                # leave it empty, but the best guess is that it is the
                # canonical transcript, accesible from 'variant' object
                # transcript_id = variant.ensembl_transcript_stable_id
                transcript_id = ''
            """
            transcript_id = gene_transcript_id

            try:
                by_transcript[transcript_id].append(variant)
            except KeyError:
                by_transcript[transcript_id] = [variant]

        variants_by_gene_by_transcript[gene] = by_transcript

    o.print('All variants', all_variants_count)

    return variants_by_gene_by_transcript, cosmic_genes_to_load


def get_all_variants(variants_by_transcript):

    unique_variants = []
    unique_variants_keys = []

    for gene_transcript_id, variants in variants_by_transcript.iteritems():
        for variant in variants:
            key = '\t'.join([
                variant.chr_name,
                variant.chrom_start,
                variant.chrom_end,
                variant.vcf_data.REF,
                ''.join([str(n) for n in variant.vcf_data.ALT]),
                variant.ensembl_gene_stable_id
            ])
            if key not in unique_variants_keys:
                variant.affected_transcripts = [gene_transcript_id]
                unique_variants.append(variant)
                unique_variants_keys.append(key)
            else:
                variant.affected_transcripts += [gene_transcript_id]

    return unique_variants


def summarize(variants_by_gene_by_transcript, cna, transcripts_positions, cosmic_mappings):

    cna_count = {}
    for gene, variants_by_transcript in variants_by_gene_by_transcript.iteritems():
        cna_count[gene] = {}
        for gene_transcript_id, variants in variants_by_transcript.iteritems():

            cna_list = cna.get_by_gene_and_transcript(gene_transcript_id)

            print(gene_transcript_id)

            transcript = cosmic_mappings.transcript_by_gene(gene_transcript_id)
            # współrzędne technicznie dobre ale inne niż cosmicowe :(
            start, end = map(int, transcripts_positions[transcript])

            gain = CopyNumberMap(start, end)
            loss = CopyNumberMap(start, end)

            for is_gain, s, e in cna_list:
                if is_gain:
                    gain.add(s, e)
                else:
                    loss.add(s, e)

            gain.show()
            loss.show()

            cna_longest_gain = gain.get_longest()
            cna_max_loss = loss.get_max()

            print('Gain/loss', cna_longest_gain, cna_max_loss)

            cna_count[gene][gene_transcript_id] = (cna_longest_gain, cna_max_loss)

            print(cna_list)

    genes_poly_aaa_increase = 0
    genes_poly_aaa_decrease = 0
    genes_poly_aaa_constant = 0
    gain_in_increase = 0
    loss_in_increase = 0
    gain_in_decrease = 0
    loss_in_decrease = 0
    loss_in_constant = 0
    gain_in_constant = 0
    variants_count = 0

    for gene, variants_by_transcript in variants_by_gene_by_transcript.iteritems():

        # opcja pierwsza - traktuj wysztkie tak samo ale wyrzuć duplikaty
        variants = get_all_variants(variants_by_transcript)

        # opcja druga - weź tylko grupę która odpowiada najdłuższemu transkryptowi:
        # variants = variant_group[] TODO TBD
        # opca trzecia - weź tylko grupę która odpowiada kanonicznemu transkryptowi
        # variants = variant_group['']

        variants_count += len(variants)

        o.print(gene)
        o.indent()
        o.print('Total:', len(variants))

        poly_a_variants = filter(lambda variant: variant.has_poly_a, variants)
        poly_a_potential_variants = filter(lambda variant: variant.will_have_poly_a, variants)

        #poly_a_constant = filter(lambda variant: not (variant.poly_a_increase or variant.poly_a_decrease), variants)
        poly_a_increase = filter(lambda variant: variant.poly_a_increase, variants)
        poly_a_decrease = filter(lambda variant: variant.poly_a_decrease, variants)

        o.print('# of poly_a adjacent to variants:')
        o.indent()
        o.print('before mutations', len(poly_a_variants))
        o.print('after mutations', len(poly_a_potential_variants))
        o.outdent()
        o.print('# of poly_a increase events', len(poly_a_increase))
        o.print('# of poly_a decrease events', len(poly_a_decrease))

        if len(poly_a_increase) > len(poly_a_decrease):
            genes_poly_aaa_increase += 1
            gain_in_increase += sum([x[0] for x in cna_count[gene].values()])
            loss_in_increase += sum([x[1] for x in cna_count[gene].values()])
        elif len(poly_a_increase) < len(poly_a_decrease):
            genes_poly_aaa_decrease += 1
            gain_in_decrease += sum([x[0] for x in cna_count[gene].values()])
            loss_in_decrease += sum([x[1] for x in cna_count[gene].values()])
        else:
            genes_poly_aaa_constant += 1
            gain_in_constant += sum([x[0] for x in cna_count[gene].values()])
            loss_in_constant += sum([x[1] for x in cna_count[gene].values()])

        o.outdent()

    gene_count = len(variants_by_gene_by_transcript)

    print('Analyzed genes:', gene_count)
    print('Analyzed variants:', variants_count)

    print('\nSummary for all genes:')
    print('\t\tCount\tCNV gain\tCNV loss')
    print('Decreased poly_aaa:', genes_poly_aaa_decrease, gain_in_decrease, loss_in_decrease)
    print('Increased poly_aaa:', genes_poly_aaa_increase, gain_in_increase, loss_in_increase)
    print('Constant poly_aaaa:', genes_poly_aaa_constant, gain_in_constant, loss_in_constant)


if __name__ == '__main__':

    import argparse

    to_show = ['databases', 'datasets', 'filters', 'attributes', 'attributes_by_page', 'some_variant']
    cache_actions = ['load', 'save']
    parser = argparse.ArgumentParser(description='Find SNPs')
    parser.add_argument('--show', choices=to_show, default=None)
    parser.add_argument('--cache', choices=cache_actions, default=None)

    subcommands = ['show', 'cache']
    arguments = ['--' + a if a in subcommands else a for a in sys.argv[1:]]

    args = parser.parse_args(arguments)

    if args.show:
        what = args.show
        if what == 'databases':
            from biomart import BiomartServer
            server = BiomartServer(BIOMART_URL)
            server.show_databases()
        if what == 'datasets':
            from biomart import BiomartServer
            server = BiomartServer(BIOMART_URL)
            server.show_datasets()
        if what == 'filters':
            som_snp = BiomartDataset(BIOMART_URL, name=DATASET)
            som_snp.show_filters()
        if what == 'attributes':
            som_snp = BiomartDataset(BIOMART_URL, name=DATASET)
            som_snp.show_attributes()
        if what == 'attributes_by_page':
            som_snp = BiomartDataset(BIOMART_URL, name=DATASET)
            som_snp.show_attributes_by_page()
        if what == 'some_variant':
            variants_by_gene, variants_data = parse(1)
            gene, variants = variants_by_gene.popitem()
            print(variants[0])
    else:
        import cPickle as pickle

        cache = args.cache
        cache_name = '.cache'

        if cache == 'load':
            with open(cache_name, 'rb') as f:
                variants_by_gene, cds_db, cdna_db = pickle.load(f)
            o.print('Variants data loaded from cache')
        else:
            o.print('Parsing dataset:', DATASET)
            variants_by_gene, variants_data = parse(100)
            transcripts_to_load = []
            for gene, variants in variants_by_gene.iteritems():
                transcripts_to_load += [v.ensembl_transcript_stable_id for v in variants]
            transcripts_to_load = set(transcripts_to_load)
            cds_db = SequenceDB(index_by='transcript', sequence_type='cds', restrict_to=transcripts_to_load)
            cdna_db = SequenceDB(index_by='transcript', sequence_type='cdna', restrict_to=transcripts_to_load)
            if cache == 'save':
                # a little trick - the iterator is a generator which could
                # not be dumped, but it is not necessary to store it at all
                it = variants_data.iterator
                variants_data.iterator = None
                with open(cache_name, 'wb') as f:
                    pickle.dump((variants_by_gene, cds_db, cdna_db), f, protocol=pickle.HIGHEST_PROTOCOL)
                o.print('variants data saved to cache')
                variants_data.iterator = it

        variants_by_gene_by_transcript, cosmic_genes_to_load = parse_variants(cds_db, cdna_db)
        from cna_by_transcript import CompleteCNA, CosmicMappings

        if cache == 'load':
            with open(cache_name + '-cna', 'rb') as f:
                cna = pickle.load(f)
            o.print('CNA loaded from cache')
        else:
            cna = CompleteCNA('cosmic/CosmicCompleteCNA.tsv', restrict_to=cosmic_genes_to_load)
            if cache == 'save':
                with open(cache_name + '-cna', 'wb') as f:
                    pickle.dump(cna, f, protocol=pickle.HIGHEST_PROTOCOL)

        cosmic_mappings = CosmicMappings()
        ensembl_transcripts = [cosmic_mappings.transcript_by_gene(gene) for gene in cosmic_genes_to_load]

        transcripts = BiomartData(
            dataset=BiomartDataset(
                BIOMART_URL,
                name='hsapiens_gene_ensembl'
            ),
            attributes=[
                'transcript_flank',
                'ensembl_transcript_id',
                'strand',
                '3_utr_start',
                '3_utr_end',
                '5_utr_start',
                '5_utr_end'
            ],
            filters={
                'upstream_flank': '1',
                'ensembl_transcript_id': ensembl_transcripts
            },
            fasta=True
        )

        transcripts_positions = {}
        for transcript in transcripts:
            header_data = transcript.header.split('|')

            transcript_id, strand = header_data[:2]
            o.print(transcript.header)
            o.indent()
            coordinates = [map(int, c.split(';')) for c in header_data[2:]]
            o.print('Coords:', coordinates)
            start_3, end_3, start_5, end_5 = coordinates

            if strand == '-1':
                end = min(start_5)
                start = max(end_3)
            else:
                end = min(start_3)
                start = max(end_5)
                assert strand == '1' or strand == '+1'

            pos = (start + 1, end - 1)

            transcripts_positions[transcript_id] = pos
            o.outdent()
            o.print('Position:', pos)

        summarize(variants_by_gene_by_transcript, cna, transcripts_positions, cosmic_mappings)
