
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

