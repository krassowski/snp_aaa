import json
from collections import defaultdict

from tqdm import trange

from commands import SourceSubparser
from cache import cacheable
from variant import Variant, BiomartVariant
from variant_sources import variants_getter

from pandas import read_csv
import numpy as np
"""
./snp_parser.py -n -d clinvar_variants clinvar --trait 'Polycystic Kidney Disease'

chrom   pos     ref     alt     measureset_type measureset_id   rcv     allele_id
symbol  hgvs_c  hgvs_p  molecular_consequence   clinical_significance
pathogenic      benign  conflicted      review_status   gold_stars
all_submitters  all_traits      all_pmids       inheritance_modes
age_of_onset     prevalence      disease_mechanism       origin  xrefs
"""

clinvar_args = SourceSubparser(
    'clinvar',
    help='Arguments for clinvar variants source'
)

clinvar_args.add_command(
    '--path',
    help='Path to clinvar tab-delimeted file',
    default='clinvar/output/b37/single/clinvar_alleles.single.b37.tsv.gz'
)

clinvar_args.add_command(
    '--trait',
    help='Restrict to variants with given trait'
)


@variants_getter
def clinvar_variants(args):

    # assert args.clinvar

    types = {
        'chrom': str,
        'pos': np.int32,
        'ref': str,
        'alt': str,
        'symbol': str,
        'hgvs_c': str,
        'hgvs_p': str,
        'molecular_consequence': str,
        'all_traits': str
    }

    df = read_csv(args.path, sep='\t', usecols=types.keys(), dtype=types)

    if args.trait:
        df = df[df.all_traits.str.contains(args.trait)]

    variants_by_gene = defaultdict(list)

    for row in df.iterrows():
        row = row[1]

        gene = row.symbol

        v = Variant(
            chr_name=row.chrom,
            chrom_start=row.pos,
            chrom_end=row.pos+len(row.alt)-len(row.ref),
            ref=row.pos,
            refsnp_id=row.hgvs_p
        )
        variants_by_gene[gene].append(v)

    return variants_by_gene
