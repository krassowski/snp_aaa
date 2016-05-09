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
import sys, os
better_biomart_path = os.path.realpath(os.path.join(os.curdir, 'biomart'))
sys.path.insert(0, better_biomart_path)

from poly_a import has_poly_a
from fasta_sequence_db import SequenceDB, FastSequenceDB
from output_formatter import OutputFormatter
from data_store import DataStore


o = OutputFormatter()

from biomart import BiomartDataset
import tabix


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


class CosmicMappings(object):

    def __init__(self, path='cosmic/All_COSMIC_Genes.fasta'):
        self.transcript_to_gene = {}
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    gene, transcript = line[1:-1].split()
                    self.transcript_to_gene[transcript] = gene

    def gene_by_transcript(self, transcript):
        return self.transcript_to_gene[transcript]


class CompleteCNA(DataStore):

    """
    a może indeksować tabixem po gene_id?
    """

    def __init__(self, path):
        header = 'ID	ID_GENE	GENE_NAME	ID_SAMPLE	ID_tumour	Primary site	Site subtype 1	Site subtype 2	Site subtype 3	Primary histology	Histology subtype 1	Histology subtype 2	Histology subtype 3	SAMPLE_NAME	TOTAL_CN	MINOR_ALLELE	MUT_TYPE	ID_STUDY	GRCh	Chromosome	G_Start	G_Stop'
        # TODO: use tabix to retrive header
        self.attributes = header.lower().split('\t')
        self.tb = tabix.open(path)
        self.cosmic = CosmicMappings()

    def get(self, chrom, start, end, transcript_id=None, flanks=0):
        start -= flanks
        end += flanks
        records = self.tb.query(chrom, start, end)

        records = map(self.parse_row, records)
        if transcript_id:
            gene_name = self.cosmic.gene_by_transcript(transcript_id)
            records = filter(lambda x: x.gene_name == gene_name, records)

        return records


class BiomartData(DataStore):

    def __init__(self, attributes, filters, dataset=None):
        self.attributes = attributes
        self.filters = filters
        if not dataset:
            dataset = BiomartDataset(BIOMART_URL, name=DATASET)
        self.dataset = dataset
        self.fetch_data()

    def fetch_data(self):
        self.data = self.dataset.search({
            'filters': self.filters,
            'attributes': [unicode(x) for x in self.attributes]
            })

        def count(self):
            return self.dataset.count()


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

        """
        Available options:

        3_prime_UTR_variant
        splice_acceptor_variant
        intergenic_variant
        UTR_variant
        coding_transcript_variant
        nonsynonymous_variant
        inframe_deletion
        downstream_gene_variant
        regulatory_region_amplification
        sequence_comparison
        TF_binding_site_variant
        transcript_amplification
        5_prime_UTR_variant
        inframe_indel
        transcript_ablation
        synonymous_variant
        sequence_variant
        terminator_codon_variant
        feature_variant
        splice_site_variant
        stop_retained_variant
        structural_variant
        protein_altering_variant
        splice_donor_variant
        inframe_insertion
        stop_lost
        feature_truncation
        feature_amplification
        exon_variant
        inframe_variant
        NMD_transcript_variant
        non_coding_transcript_exon_variant
        splice_region_variant
        TFBS_ablation
        transcript_variant
        incomplete_terminal_codon_variant
        stop_gained
        coding_sequence_variant
        gene_variant
        upstream_gene_variant
        regulatory_region_ablation
        TFBS_amplification
        start_lost
        frameshift_variant
        regulatory_region_variant
        splicing_variant
        feature_elongation
        missense_variant
        mature_miRNA_variant
        intron_variant
        internal_feature_elongation
        feature_ablation
        non_coding_transcript_variant
        """
        filters.update({
            'so_parent_name':
                [
                    'synonymous_variant',
                    'missense_variant',
                    'stop_gained',
                    'coding_sequence_variant'
                ]
        })
        super(self.__class__, self).__init__(attributes, filters, dataset)


def show_pos_with_context(seq, start, end):
    return seq[:start] + '→' + seq[start:end] + '←' + seq[end:]


def gene_names_from_patacsdb_csv(how_many=250):
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


def parse():

    genes_from_patacsdb = gene_names_from_patacsdb_csv()

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
    for line in variants.data.iter_lines():
        variant = variants.parse_row(line)
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

    return variants_by_gene


def ref_seq_len(src, ref):
    if src not in ref:
        return 0
    return len(ref[src].strip('-'))


def analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl, cna):
    o.mute()

    o.print('Variant name:', variant.refsnp_id)
    o.indent()

    # gene_id = variant.ensembl_gene_stable_id

    # trasnscript is a transcript of reference, not of the variant!
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
    # TODO: remove plus one in fetch and showpostwith context call

    pos = [str(variant.chr_name), int(variant.chrom_start), int(variant.chrom_end)]

    seq = chromosome.fetch(pos[1], pos[2], offset)
    # nuc = chromosome.fetch(pos[1], pos[2])

    reference_nuc['genome'] = seq[offset:-offset]
    reference_seq['genome'] = seq

    o.print('Chromosomal position: {0}:{1}-{2}'.format(*pos))


    # to powinno być trzymane na poziomie transkryptu?
    #   cna_records = cna.get(*pos, transcript_id=variant.ensembl_transcript_stable_id, flanks=10)

    # skoro one są takie duże to nie opłaca mi się tego indexowac na współrzędnych
    # mogę spróbować agregować to przy wczytywaniu → najpierw mam warianty, potem wczytuję:
    # biorę nazwę genu, zapisuję koordynaty oraz gain/loss. copy number już mnie nie interesuje.
    # robię algorytm grupujący gain/loss i względem koordynatów, tak jak na stronie. potem wybieram
    # odpowiednio: # of najdłuższy gain, max(loss) i zapisuję.
    # robię to niezależnie dla każdego "genu cosmicowego", potem grupuję to w prawdziwe "geny". i patrzę czy konsystentne
    # - wzgledem siebie i jako całość - szczególnie stosunki gain/loss muszą być zbliżone!

    # tak samo z rekordami z ensembla - cosmic grupuję nie dość że po genie prawdziwym to po genie "cosmicowym" i sprawdzam czy ilość
    # mutacji poliA jest konsystentna. Jeśli nie uda się tego sensownie zrobić, to po prostu je grupuję.

    #   o.print("CNA records count: {}".format(len(cna_records)))

    #   variant.cna_list = cna_records

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
        o.print('Lack of variant data! Abort!')
        o.print(variant.refsnp_id)
        # eg:
        # cancer.sanger.ac.uk/cosmic/mutation/overview?id=250006
        # cancer.sanger.ac.uk/cosmic/mutation/overview?id=4967608
        variant.correct = False
        # ustalenie: pomijać

        o.outdent()
        return False

    #o.unmute()
    # o.print(record.num_called, record.call_rate, record.num_unknown)  # ZeroDivisionError ?
    # o.print(record.num_hom_ref, record.num_het, record.num_hom_alt)
    # o.print(record.nucl_diversity, record.aaf, record.heterozygosity)  # ZeroDivisionError ?
    # o.print(record.is_snp, record.is_indel, record.is_transition, record.is_deletion)
    # o.print(record.is_monomorphic)
    # o.print(vcf_data)
    #exit()

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


def summarize(variants_by_gene, cds_db, cdna_db, cna):

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

    for gene, variants in variants_by_gene.iteritems():

        # Just to be certain
        variants_unique_ids = set(variant.refsnp_id for variant in variants)
        assert len(variants_unique_ids) == len(variants)

        all_variants_count += len(variants)

        for variant in variants:
            analyze_variant(variant, cds_db, cdna_db, dna_db, vcf_cosmic, vcf_ensembl, cna)

        # Remove variants with non-complete data 
        correct_variants = filter(lambda variant: variant.correct, variants)

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

            full_gene_id = variant.vcf_data.INFO['GENE']
            chunks = full_gene_id.split('_')
            if len(chunks) > 1:
                transcript_id = chunks[1]
                # just in case to detect potential ambuguity
                assert len(chunks) == 2
            else:
                # leave it empty, but the best guess is that it is the
                # canonical transcript, accesible from 'variant' object
                # transcript_id = variant.ensembl_transcript_stable_id
                transcript_id = ''

            try:
                by_transcript[transcript_id].append(variant)
            except KeyError:
                by_transcript[transcript_id] = [variant]

        variants_by_gene_by_transcript[gene] = by_transcript

    genes_poly_aaa_increase = 0
    genes_poly_aaa_decrease = 0
    copy_number_in_increase = 0
    copy_number_in_decrease = 0
    variants_count = 0

    for gene, variants_by_transcript in variants_by_gene_by_transcript.iteritems():

        # opcja pierwsza - traktuj wysztkie tak samo ale wyrzuć duplikaty
        unique_variants = []
        unique_variants_keys = []

        for transcript_id, variants in variants_by_transcript.iteritems():
            # Loss of information here!
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
                    unique_variants.append(variant)
                    unique_variants_keys.append(key)

        variants = unique_variants

        # variants = linearize(variant_groups)
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

        poly_a_increase = filter(lambda variant: variant.poly_a_increase, variants)
        poly_a_decrease = filter(lambda variant: variant.poly_a_decrease, variants)

        o.print('# of poly_a before mutations', len(poly_a_variants))
        o.print('# of poly_a after mutations', len(poly_a_potential_variants))
        o.print('# of poly_a increase events', len(poly_a_increase))
        o.print('# of poly_a decrease events', len(poly_a_decrease))

        if len(poly_a_increase) > len(poly_a_decrease):
            genes_poly_aaa_increase += 1
        elif len(poly_a_increase) < len(poly_a_decrease):
            genes_poly_aaa_decrease += 1

        o.outdent()

    gene_count = len(variants_by_gene_by_transcript)

    print('Analyzed genes:', gene_count)
    print('Analyzed variants:', variants_count)
    print('Rejected variants:', all_variants_count - variants_count)
    print('Decreased poly_aaa:', genes_poly_aaa_decrease)
    print('Increased poly_aaa:', genes_poly_aaa_increase)

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
            o.print('Parsing dataset:', DATASET)
            variants_by_gene = parse()
            transcripts_to_load = []
            for gene, variants in variants_by_gene.iteritems():
                transcripts_to_load += [v.ensembl_transcript_stable_id for v in variants]
            transcripts_to_load = set(transcripts_to_load)
            cds_db = SequenceDB(index_by='transcript', sequence_type='cds', restrict_to=transcripts_to_load)
            cdna_db = SequenceDB(index_by='transcript', sequence_type='cdna', restrict_to=transcripts_to_load)
            if cache == 'save':
                with open(cache_name, 'wb') as f:
                    pickle.dump((variants_by_gene, cds_db, cdna_db), f, protocol=pickle.HIGHEST_PROTOCOL)
                o.print('variants data saved to cache')

        # tego będę używać
        #cna = CompleteCNA2('cosmic/CosmicCompleteCNA.tsv', restrict_to=transcripts_to_load)
        # to tylko chwilowe:
        cna = CompleteCNA('cosmic/CosmicCompleteCNA.tsv.sorted.gz')

        summarize(variants_by_gene, cds_db, cdna_db, cna)

