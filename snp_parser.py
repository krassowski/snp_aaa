#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

BIOMART_URL = 'http://www.ensembl.org/biomart'
"""
Mirrors:
    uswest.ensembl.org
    useast.ensembl.org
    asia.ensembl.org
"""
DATASET = 'hsapiens_snp_som'

#pip install pyvcf pysam biomart
import vcf
# pysam is required too
import sys
#better_biomart_path = os.path.realpath(os.path.join(os.curdir, 'biomart'))
#sys.path.insert(0, better_biomart_path)

from poly_a import has_poly_a
from fasta_sequence_db import SequenceDB, FastSequenceDB
from output_formatter import OutputFormatter

o = OutputFormatter()

from biomart import BiomartDataset


class BiomartData(object):

    def __init__(self, attributes, filters):
        self.attributes = attributes
        self.filters = filters
        self.fetch_data()

    def fetch_data(self):
        self.data = som_snp.search({
            'filters': self.filters,
            'attributes': [unicode(x) for x in self.attributes]
        })

    def count(self):
        return som_snp.count()

    def pos(self, attribute):
        try:
            return self.attributes.index(attribute)
        except ValueError:
            raise KeyError

    def parse_row(self, row):
        return DataRow(self, row)


class DataRow(object):

    def __init__(self, parent, line):
        line = line.decode('utf-8')
        self._data = line.split('\t')
        self._parent = parent

    def __getattr__(self, k):
        if k.startswith('__') and k.endswith('__'):
            return super(DataRow, self).__getattr__(k)
        try:
            return self._data[self._parent.pos(k)]
        except KeyError:
            raise AttributeError

    def __setattr__(self, k, v):
        if k.startswith('_'):
            super(DataRow, self).__setattr__(k, v)
        else:
            try:
                self._data[self.get_index(k)] = v
            except KeyError:
                self._parent.attributes.append(k)
                self.__setattr__(k, v)
            except IndexError:
                while len(self._data) <= len(self._parent.attributes):
                    self._data.append(None)
                self.__setattr__(k, v)

    def __repr__(self):
        buff = 'Biomart DataRow:\n\t'
        buff += '\n\t'.join([self.get_key(i) + ': ' + str(v) for i, v in enumerate(self._data)])

        return buff

    def get_key(self, index):
        return self._parent.attributes[index]

    def get_index(self, key):
        return self._parent.pos(key)

    def get(self, k):
        return self.__getattr__(k)

    def set(self, k, v):
        return self.__setattr__(k, v)


class VariantsData(BiomartData):

    def __init__(self, attributes=None, filters=None):

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

        # TODO być może to właśnie ogranicza do cosmica
        # synonymous_variant
        # lista variantów na mailu
        # 'so_parent_name': ['coding_sequence_variant', 'synonymous_variant']
        # problemy z biomartem..
        filters.update({
            #'so_parent_name': 'exon_variant'
            'so_parent_name': 'coding_sequence_variant'
        })
        super(self.__class__, self).__init__(attributes, filters)


def show_pos_with_context(seq, start, end):
    return seq[:start] + '→' + seq[start:end] + '←' + seq[end:]


def gene_names_from_patacsdb_tsv():
    genes = []
    with open('patacsdb_first_100.tsv', 'r') as f:
        for line in f:
            if line:
                data = [x.strip() for x in line.split('\t')]
                genes.append(data[1])
    return genes


def parse():

    genes_from_patacsdb = gene_names_from_patacsdb_tsv()

    filters = {
        'ensembl_gene': genes_from_patacsdb
    }

    variants = VariantsData(filters=filters)

    variants_by_gene = {}

    previous_id = None
    # te atrybuty są zmienną częścią w zduplikowanych rekordach -
    # dla jednej mutacji o danym id, ensembl zwraca np 4 rekordy
    # z których w każdym wszystkio jest identyczne oprócz poniższych atrybutów
    # (testowane na cosmicach) - tak więc zczytując kompresuję to do jednego wpisu
    # w którym trzymam te wartości w postaci uporzadkowanej listy
    attrs = ['phenotype_description', 'set_name', 'set_description']
    for line in variants.data.iter_lines():
        variant = variants.parse_row(line)
        gene = variant.ensembl_gene_stable_id
        if variant.refsnp_id == previous_id:
            last = variants_by_gene[gene][-1]
            assert last.refsnp_id == previous_id
            for attr in attrs:
                last.get(attr).append(variant.get(attr))
        else:
            for attr in attrs:
                variant.set(attr, [variant.get(attr)])
            previous_id = variant.refsnp_id
            if gene not in variants_by_gene:
                variants_by_gene[gene] = []
            variants_by_gene[gene].append(variant)

    return variants_by_gene

def ref_seq_len(src, ref):
    if not src in ref:
        return 0
    return len(ref[src].strip('-'))


def analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl):
    o.mute()

    o.print('Variant name: ' + variant.refsnp_id)
    o.indent()

    gene_id = variant.ensembl_gene_stable_id

    transcript_id = variant.ensembl_transcript_stable_id
    strand = int(variant.ensembl_transcript_chrom_strand)

    offset = 20

    o.print('Transcript ' + transcript_id)
    o.indent()

    if variant.refsnp_source != 'COSMIC':
        o.unmute()
    o.print(variant.consequence_type_tv)
    o.print('strand: ' + str(strand))

    reference_nuc = {}
    reference_seq = {}

    for db in [cdna_db, cds_db]:
        src = db.sequence_type

        start, end = variant.get(src + '_start'), variant.get(src + '_end')
        seq = db.fetch(transcript_id, strand, start, end, offset)

        if not db.has(transcript_id) or not seq:
            o.print('Lack of transcript in ' + src)
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
    # todo: remove plus one in fetch and showpostwith context call

    pos = [str(variant.chr_name), int(variant.chrom_start), int(variant.chrom_end)]

    seq = chromosome.fetch(pos[1], pos[2], offset)
    #nuc = chromosome.fetch(pos[1], pos[2])

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
    for src, sequence in reference_seq.items():
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
        # todo: skip?
        o.print('Reference sequences are not consistent:')
        for src, seq in reference_seq.items():
            if seq:
                o.print(src + ':\t' + show_pos_with_context(seq, offset, -offset))

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

    if ref_seq_len('cdna', reference_seq) > ref_seq_len('cds', reference_seq) and ref_seq_len('cdna', reference_seq):
        ref_seq = reference_seq['cdna']
        o.print('\tChosing cdna sequence as reference')
    elif ref_seq_len('cds', reference_seq) > ref_seq_len('cdna', reference_seq) and ref_seq_len('cds', reference_seq):
        ref_seq = reference_seq['cds']
        o.print('\tChosing cds sequence as reference')
    else:
        ref_seq = reference_seq['genome']
        o.print('\tChosing genome sequence as reference')


    o.print('Context: ' + show_pos_with_context(ref_seq, offset, -offset))

    main_vcf = vcf_cosmic if variant.refsnp_source == 'COSMIC' else vcf_ensembl

    vcf_data = None

    for record in main_vcf.fetch(*pos):

        if record.ID == variant.refsnp_id:
            vcf_data = record
            break

    if not vcf_data:
        o.unmute()
        o.print('Lack of varaiant data! Abort!')
        o.print(variant.refsnp_id)
        # eg:
        # cancer.sanger.ac.uk/cosmic/mutation/overview?id=250006
        # cancer.sanger.ac.uk/cosmic/mutation/overview?id=4967608
        variant.correct = False
        # ustalenie: pomijać

        o.outdent()
        return False

    # o.print(record.num_called, record.call_rate, record.num_unknown)  # ZeroDivisionError ?
    # o.print(record.num_hom_ref, record.num_het, record.num_hom_alt)
    # o.print(record.nucl_diversity, record.aaf, record.heterozygosity)  # ZeroDivisionError ?
    # o.print(record.is_snp, record.is_indel, record.is_transition, record.is_deletion)
    # o.print(record.is_monomorphic)
    # o.print(vcf_data)

    alt = vcf_data.ALT
    if len(alt) != 1:
        o.print('Too much alternative alleles')

    mutated_seq = ref_seq[:offset] + str(alt[0]) + ref_seq[offset + 1:]

    variant.sequence = ref_seq
    variant.vcf_data = vcf_data
    variant.has_poly_a = has_poly_a(ref_seq, offset, len(ref_seq) - offset)
    variant.will_have_poly_a = has_poly_a(mutated_seq, offset, len(mutated_seq) - offset)
    variant.correct = True

    o.print(mutated_seq)

    o.unmute()
    o.outdent()

    if variant.allele != 'COSMIC_MUTATION':
        exit()

    return True

def summarize(variants_by_gene, cds_db, cdna_db):

    """
    vcf.parser uses 0-based coordinates:
    http://pyvcf.readthedocs.org/en/latest/_modules/vcf/parser.html?highlight=coordinates
    ensembl uses 1-based coordinates:
    http://www.ensembl.org/info/docs/api/core/core_tutorial.html#coordinates
    """

    chromosomes = map(str, range(1, 22)) + ['X', 'Y']

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

    for gene, variants in variants_by_gene.iteritems():
        o.print(gene)
        o.indent()
        o.print('Total: ' + str(len(variants)))

        unique_variants = set(variant.refsnp_id for variant in variants)
        assert len(unique_variants) == len(variants)

        for variant in variants:
            analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl)

        correct_variants = filter(lambda variant: variant.correct, variants)
        poly_a_variants = filter(lambda variant: variant.has_poly_a, variants)
        poly_a_potentiall_variants = filter(lambda variant: variant.will_have_poly_a, variants)

        o.print(len(poly_a_variants))
        o.print(len(poly_a_potentiall_variants))

        o.outdent()

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
            variants_by_gene = parse()
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
            o.print('Parsing dataset: ' + DATASET)
            som_snp = BiomartDataset(BIOMART_URL, name=DATASET)
            variants_by_gene = parse()
            variants_to_load = []
            for gene, variants in variants_by_gene.iteritems():
                variants_to_load += [v.ensembl_transcript_stable_id for v in variants]
            variants_to_load = set(variants_to_load)
            cds_db = SequenceDB(index_by='transcript', sequence_type='cds', restrict_to=variants_to_load)
            cdna_db = SequenceDB(index_by='transcript', sequence_type='cdna', restrict_to=variants_to_load)
            if cache == 'save':
                with open(cache_name, 'wb') as f:
                    pickle.dump((variants_by_gene, cds_db, cdna_db), f)
                o.print('variants data saved to cache')


        summarize(variants_by_gene, cds_db, cdna_db)

