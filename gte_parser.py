#!/usr/bin/python3
import glob
from multiprocessing import Pool


def get_files(path, pattern):
    return glob.glob(path + '/' + pattern)


def get_gbff_files(path='./ncbi/mRNA_Prot', pattern='human.*.rna.gbff'):
    return get_files(path, pattern)


def get_gtex_files(path='./gtex/GTEx_Analysis_V6_eQTLs', pattern='*.snpgenes'):
    return get_files(path, pattern)

"""
Interesują mnie tylko takie rekordy w gbff które są kodujące:
        mol_type="mRNA"
    równoważnie
        translation=
"""


class Gene:

    __slots__ = 'accession', 'seq', 'variants'

    def __init__(self):
        self.seq = ''
        self.variants = []


class Variant:

    __slots__ = 'pos', 'ref', 'alt', 'rs', 'compl', 'len'

    def __init__(self, pos, rs, compl, length, ref, alt):
        self.pos = pos
        self.ref = ref
        self.compl = compl
        self.alt = alt
        self.rs = rs
        self.len = length


def get_genes_from_file(filename):
    """
    Pomysł: weź accession, potem je wszystkie wyślesz do biomarta i
    przekonwertujesz na ensemblowe
    https://biodbnet-abcc.ncifcrf.gov/db/db2dbRes.php
    """
    genes = []

    gene = Gene()
    accepted = False

    with open(filename, 'r') as f:

        for line in f:

            line = line.strip()

            if line.startswith('ACCESSION'):
                gene.accession = line.split()[1]
                continue

            if line.startswith('variation'):
                length = 1

                pos = line.split()[-1]

                compl = False

                if pos.startswith('complement'):
                    pos = pos[11:-1]
                    compl = True

                if '..' in pos:
                    d = tuple(map(int, pos.split('..')))
                    length  = d[1] - d[0] + 1
                    pos = d[0]



                pos = int(pos)

                # ref, alt
                alleles = []
                rs = None

                allel_opened = False

                while True:
                    line = next(f).strip()

                    if line.startswith('/'):
                        data = line.split('=')
                        value = data[1].strip('"')
                        if data[0] == '/replace':
                            allel_opened = line.count('"') == 1
                            if len(alleles) < 2:
                                alleles.append(value)
                            #else:
                            #    print('More than two alleles!')

                        if data[0] == '/db_xref':
                            assert value.startswith('dbSNP:')
                            rs = value[6:]
                            break
                    else:
                        #print(gene.accession)
                        assert rs is None
                        if allel_opened:
                            allel_opened = line.count('"') == 1
                            alleles[-1] += line.strip('"')

                try:
                    assert len(alleles) == 2
                except:
                    print (alleles, rs)

                variant = Variant(pos, rs, compl, length, *alleles)

                gene.variants.append(variant)
                continue

            if line.startswith('source'):

                next(f)    # /organism
                line = next(f).strip()

                assert '/mol_type=' in line

                mol_type = line.split('=')[1]
                #print(mol_type)
                if mol_type == '"mRNA"':
                    accepted = True
                    continue

            if line.startswith('ORIGIN'):

                line = next(f)
                while line != '//\n':
                    line = line[10:]
                    gene.seq += ''.join(line.split())
                    line = next(f)

                if accepted:
                    genes.append(gene)
                    gene = Gene()
                continue

    return genes


def check_variants(gene):
    for variant in gene.variants:
        seq = gene.seq
        if variant.compl:
            seq = complement(seq)[::-1]
        print (seq[variant.pos - 1:variant.pos -1 + variant.len], variant.ref)
        assert seq[variant.pos - 1:variant.pos -1 + variant.len] == variant.ref


def gather_coding_variations_per_gene():
    parsing_pool = Pool()
    filenames = get_gbff_files()[:2]
    genes = parsing_pool.map(get_genes_from_file, filenames)
    return genes


def complement(seq):
    basic = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    try:
        return ''.join([basic[n] for n in seq])
    except KeyError:
        # http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
        IUPAC = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A', 'Y': 'R',
                 'R': 'Y', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
        return ''.join([IUPAC[n] for n in seq])


def get_all_cds_variants(attributes=None, filters=None):

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
        'minor_allele',
        'chrom_strand',
        'cdna_start',
        'cdna_end',
        'ensembl_gene_stable_id',
        'ensembl_transcript_stable_id',
        'ensembl_transcript_chrom_strand',
        'cds_start',
        'cds_end'
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

    return self.dataset.search({
        'filters': self.filters,
        'attributes': [unicode(x) for x in self.attributes]
        })

