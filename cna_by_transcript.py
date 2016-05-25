# -*- coding: utf-8 -*-


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
    # join -1 2 -2 14 <(cat pnma1.temp.tsv | sort -k2) \
    # <(cat cna.pnma1.temp.tsv | sort -k14)

    expr_map = {
        'over': 0,
        'normal': 1,
        'under': 2
    }

    def __init__(self, path, restrict_to=None, include_all=False):
        """
        include_all: if False - use 'high value (numeric)', else use 'all' mode
        """
        """
        Header:
        ID  ID_GENE GENE_NAME   ID_SAMPLE   ID_tumour   Primary site
        Site subtype 1  Site subtype 2  Site subtype 3  Primary histology
        Histology subtype 1 Histology subtype 2 Histology subtype 3 SAMPLE_NAME
        TOTAL_CN    MINOR_ALLELE    MUT_TYPE    ID_STUDY    GRCh
        Chromosome:G_Start..G_Stop
        """

        self.restrict_to = restrict_to

        db = {}

        rstrip = str.rstrip
        split = str.split

        expression_data = {}
        # parse expression file
        with open('cosmic/CosmicCompleteGeneExpression.tsv', 'r') as f:
            header = f.next().rstrip().split('\t')
            for line in f:
                line = split(rstrip(line), '\t')
                gene = line[2]

                if self.restrict_to and gene not in self.restrict_to:
                    continue

                sample = line[1]
                expr = self.expr_map[line[3]]

                expression_data[(gene, sample)] = expr

        # parse CNA file
        with open(path, 'r') as f:
            header = f.next().rstrip().split('\t')
            self.original_header = header

            for line in f:
                line = split(rstrip(line), '\t')

                # If we are in 'high value (numeric)' mod (indicated by False
                # on 'include_all') then we want to exclude lines that does
                # not have properly defined 'copy number' or 'minor allele'
                if not (include_all or (line[14] and line[15])):
                    continue

                gene = line[2]

                if self.restrict_to and gene not in self.restrict_to:
                    continue

                sample = line[13]

                try:
                    expr = expression_data[(gene, sample)]
                except KeyError:
                    # if we do not have expression data - skip the CNV
                    continue

                try:
                    db[gene][expr] += 1
                except KeyError:
                    db[gene] = [0, 0, 0]
                    db[gene][expr] += 1

        self.db = db

    def get_by_gene_and_transcript(self, gene_transcript):
        return self.db[gene_transcript]
