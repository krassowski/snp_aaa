# -*- coding: utf-8 -*-

class Chromosome(object):

    def __init__(self, name):
        self.name = name

class CosmicMappings(object):

    def __init__(self, path='cosmic/All_COSMIC_Genes.fasta'):
        self.transcript_to_gene = {}
        self.gene_to_transcript = {}
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    gene, transcript = line[1:-1].split()
                    self.transcript_to_gene[transcript] = gene
                    self.gene_to_transcript[gene] = transcript

    def gene_by_transcript(self, transcript):
        return self.transcript_to_gene[transcript]

    def transcript_by_gene(self, gene):
        return self.gene_to_transcript[gene]

class CompleteCNA(object):

    def __init__(self, path, restrict_to=None, include_all=False):
        """
        include_all: if False - use 'high value (numeric)', else use 'all' mode
        """
        """
        Header:
        ID  ID_GENE GENE_NAME   ID_SAMPLE   ID_tumour   Primary site
        Site subtype 1  Site subtype 2  Site subtype 3  Primary histology
        Histology subtype 1 Histology subtype 2 Histology subtype 3 SAMPLE_NAME
        TOTAL_CN    MINOR_ALLELE    MUT_TYPE    ID_STUDY    GRCh    Chromosome:G_Start..G_Stop
        """


        self.restrict_to = restrict_to

        db = {}
        gene_to_gene_id = {}

        rstrip = str.rstrip
        split = str.split

        chromosomes_map = map(str, range(1, 23)) + ['X', 'Y']
        chromosomes = map(Chromosome, chromosomes_map)
        chromosomes_map = {chrom: nr for nr, chrom in enumerate(chromosomes_map)}

        with open(path, 'r') as f:
            header = f.next().rstrip().split('\t')
            self.original_header = header

            # gain_loss = {'GAIN': True, 'LOSS': False}
            for line in f:
                line = rstrip(line).split('\t')

                # If we are in 'high value (numeric)' mode, which is indicated by
                # False on 'include_all', then we want to exclude lines that does
                # not have properly defined 'copy number' or 'minor allele'
                if not (include_all or (line[14] and line[15])):
                    continue

                gene_name = line[2]

                if self.restrict_to and gene_name not in self.restrict_to:
                    continue

                chrom, cords = split(line[19], ':')
                start, end = split(cords, '..')

                if line[16] == 'GAIN':
                    gain = True
                else:
                    gain = False
                    assert line[16] == 'LOSS'

                try:
                    db[gene_name][1].append((gain, int(start), int(end)))
                except KeyError:
                    db[gene_name] = (chromosomes[chromosomes_map[chrom]], [(gain, int(start), int(end))])

        self.db = db

    def get_by_gene_and_transcript(self, gene_transcript):
        return self.db[gene_transcript][1]

