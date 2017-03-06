#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import sys

biomart_fork_path = os.path.realpath(os.path.join(os.curdir, 'biomart'))
sys.path.insert(0, biomart_fork_path)

import traceback
from collections import defaultdict

from biomart import BiomartServer
from tqdm import tqdm
from tqdm import trange

from biomart_data import BiomartData
from biomart_data import BiomartDataset
from output_formatter import OutputFormatter
from variant import Variant
from cache import cacheable
from fasta_sequence_db import SequenceDB, FastSequenceDB
from commands import execute_commands
from commands import append_commands

o = OutputFormatter()

VERBOSITY_LEVEL = 0

GRCH_VERSION = 'GRCh37'
GRCH_SUBVERSION = '13'
ENSEMBL_VERSION = '75'
COSMIC_VERSION = '79'
DBSNP_VERSION = '149'
SPIDEX_LOCATION = 'spidex_public_noncommercial_v1.0/spidex_public_noncommercial_v1_0.tab.gz'

DEFAULT_BIOMART_URL = 'http://grch37.ensembl.org/biomart'


vcf_locations = {
    #'COSMIC': 'cosmic/v' + COSMIC_VERSION + '/CosmicCodingMuts.vcf.gz.bgz',
    'dbSNP': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
    GRCH_SUBVERSION + '/00-All.vcf.gz',
    'ClinVar': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
    GRCH_SUBVERSION + '/00-All.vcf.gz',
    'ensembl': 'ensembl/v' + ENSEMBL_VERSION + '/Homo_sapiens.vcf.gz'
}


class VariantsData(BiomartData):

    def __init__(self, dataset=None, filters=None):

        if filters is None:
            filters = {}

        filters.update({
            'so_parent_name':
                [
                    'protein_altering_variant',   # inframe + frameshift
                    'synonymous_variant',
                    'missense_variant',
                    'frameshift_variant'
                ]
        })
        super(self.__class__, self).__init__(dataset, Variant, filters)


def gene_names_from_patacsdb_csv(limit_to=None):
    """
    Load gene names from given csv file. The default file was exported from
    sysbio.ibb.waw.pl/patacsdb database, from homo sapiens dataset.

    The gene names are expected to be located at second column in each row.
    The count of gene names to read can be constrained with use of an optional
    `limit_to` argument (default None denotes that all gene names
    from the file should be returned)
    """
    genes = set()
    count = 0
    with open('patacsdb_all.csv', 'r') as f:
        for line in f:
            if limit_to and count >= limit_to:
                break
            count += 1
            if line:
                data = [x.strip() for x in line.split(',')]
                genes.add(data[1])

    return list(genes)


@cacheable
def get_all_human_genes(biomart):
    """Retrieves list of human gene identifiers from Ensembl's biomart."""
    gene_names = set()

    print('Downloading list of all human genes in Ensembl...')

    genes_dataset = BiomartDataset(biomart, name='hsapiens_gene_ensembl')
    response = genes_dataset.search({'attributes': ['ensembl_gene_id']})

    for line in response.iter_lines():
        gene = line.decode('utf-8')
        gene_names.add(gene)

    print('Download of human genes list finished.')

    return gene_names


def compare_variants(variant, known_variant):
    """Make sure that two variants are essentially the same,

    except for transcript ids which are allowed to differ between
    otherwise identical variants. If transcripts differ,
    transcript id from 'variant' will be appended to set of
    'affected transcripts' in 'known variant'.

    Returns:
        True if variants are identical,
        False if variants differ by transcript id only

    Raises:
        ValueError if variants differ in anything but transcript id
    """
    variable_attributes = {'ensembl_transcript_stable_id'}
    const_attributes = set(known_variant.attributes) - variable_attributes

    for attr in const_attributes:
        if getattr(known_variant, attr) != getattr(variant, attr):
            raise ValueError(
                'Received variants with the same id: %s '
                'but different %s. Only transcripts are '
                'allowed to differ in Biomart fetched data'
                %
                (variant.refsnp_id, attr)
            )

    all_the_same = True

    # this is less important check
    # (if variants are the same, it's just a duplicate - no worries)
    for attr in variable_attributes:
        if getattr(known_variant, attr) != getattr(variant, attr):
            all_the_same = False
            break

    if all_the_same:
        print(
            'Received identical variants with id: %s'
            %
            variant.refsnp_id
        )
        return True

    known_variant.affected_transcripts.update([
        variant.ensembl_transcript_stable_id
    ])

    return False


@cacheable
def download_variants(biomart, dataset, gene_names, step_size=50, filters={}):
    """Retrieve from Ensembl's biomart all variants that affect given genes.

    Variants that occur repeatedly in a given gene (i.e. were annotated for
    multiple transcripts) will be reported only once; list of affected
    transcripts will be preserved in 'affected_transcripts' property.

    The genes should be specified as the list of ensembl gene identifiers.


    Returns:
        dict of lists where variants are grouped by ensembl_gene_stable_id

        lists of variants are guaranteed to contain at least one variant
        each - if there is a key with id of given gene, it has at least one
        variant in list, fetched from Ensembl Biomart.
    """

    variants_by_gene = defaultdict(list)
    variants_count = 0

    print('Downloading variants data from Ensembl\'s Biomart:')

    if gene_names == ['all_human_genes']:
        gene_names = get_all_human_genes.load_or_create(biomart)

    for start in trange(0, len(gene_names), step_size):

        query_filters = {'ensembl_gene': gene_names[start:start + step_size]}
        query_filters.update(filters)

        print('Downloading variants for genes:', gene_names[start:start + step_size])
        variants = VariantsData(filters=filters, dataset=dataset)
        print('Download completed. Parsing...')

        variants_in_gene = 0

        for variant in variants:
            assert variant.refsnp_id
            variants_in_gene += 1

            gene = variant.ensembl_gene_stable_id

            variants_with_the_same_id = [
                known_variant
                for known_variant in variants_by_gene[gene]
                if known_variant.refsnp_id == variant.refsnp_id
            ]

            # if we have more than one, we have a duplicate in old data.
            assert len(variants_with_the_same_id) <= 1

            if variants_with_the_same_id:
                known_variant = variants_with_the_same_id[0]
                compare_variants(variant, known_variant)
            else:
                variants_by_gene[gene].append(variant)

        # this is a small trick to turn off unpickable iterator so it is
        # possible to save the variants object as a cache by pickling
        variants.iterator = None

        variants_count += variants_in_gene
        print('Parsed %s variants.' % variants_in_gene)

    print('Downloaded %s variants' % variants_count)
    return variants_by_gene


def select_poly_a_related_variants(variants):
    """Return list oof variants occurring in a poly(A) track and variants which will create poly(A) track."""
    return [
        variant
        for variant in variants
        if any([
            data.has or data.will_have
            for data in variant.poly_aaa.values()
        ])
    ]


def all_poly_a_variants(variants_by_gene):
    total = len(variants_by_gene)

    for gene, variants in tqdm(variants_by_gene.iteritems(), total=total):

        poly_a_related_variants = select_poly_a_related_variants(variants)

        for variant in poly_a_related_variants:
            yield variant


def perform_analyses(args):
    """The main workflow happens here."""

    variants_by_gene = variants_by_gene_parsed.load()
    from analyses import REPORTERS

    for reporter_name in args.report:
        current_reporter = REPORTERS[reporter_name]

        try:
            current_reporter(variants_by_gene)

        except Exception:
            traceback.print_exc()


def create_arg_parser():
    import argparse
    import json
    from analyses import REPORTERS
    from analyses import VARIANTS_GETTERS

    to_show = [
        'databases', 'datasets', 'filters', 'attributes',
        'attributes_by_page', 'some_variant'
    ]

    parser = argparse.ArgumentParser(description=(
        'Retrieve all data about variants from genes belonging to a predefined '
        'or user-chosen list of genes (see --genes_list) [including list of all '
        'human genes] and perform one or multiple of available analyses '
        '(see --report).\n'
        'Once variants are downloaded and pre-parsed, they will be stored '
        'in cache until another list of variants will be specified with '
        'options: --download_variants --parse_variants.'
    ))
    parser.add_argument(
        '--report',
        nargs='+',
        choices=REPORTERS.keys(),
        default=[],
        help='Analyses to be performed; one or more from: ' + ', '.join(REPORTERS.keys()),
        metavar=''
    )
    # parser.add_argument(
    #     '--variants_list',
    #     nargs='+',
    #     type=str,
    #     default=None,
    #     help=(
    #         'list of variants to be used in analysis. By default, '
    #         'analysis will be performed using variants fetched with biomart '
    #         'from coding areas of human genes from PATACSDB.'
    #     )
    # )
    parser.add_argument(
        '--filters',
        type=json.loads,
        default={},
        help='Additional variants filters to be used when querying biomart'
    )
    parser.add_argument(
        '--genes_list',
        nargs='+',
        type=str,
        default=gene_names_from_patacsdb_csv(),
        help=(
            'list of human genes from which variants should be extracted for use in '
            'available analyses. By default all human genes from PATACSDB '
            'will be used. To use all human genes specify "all_human_genes".'
        )
    )
    parser.add_argument(
        '--dataset',
        type=str,
        help='name of biomart dataset to be used, eg. hsapiens_snp',
        default='hsapiens_snp_som'
    )
    parser.add_argument(
        '--biomart',
        type=str,
        help=(
            'URL of biomart to be used. Default: %s. '
            'To query GRCh38 use http://www.ensembl.org/biomart/ '
            'For mirrors of GRCh38 replace "www" with: "uswest", "useast" or "asia". '
            %
            DEFAULT_BIOMART_URL
        ),
        default=DEFAULT_BIOMART_URL
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='increase output verbosity'
    )
    parser.add_argument(
        '-d', '--download_variants',
        choices=VARIANTS_GETTERS.keys(),
        default=None,
        help='How should the variants be retrieved? Choices are: '
             'biomart (downloads variants from given --biomart, '
             'belonging to coding sequences of genes specified with '
             '--genes_list and filtered with --filters), ' + ', '.join(VARIANTS_GETTERS.keys())
    )
    parser.add_argument(
        '--parse_variants',
        action='store_true',
        help='preparse variants'
    )
    parser.add_argument(
        '-s',
        '--step_size',
        type=int,
        default=25,
        help='When using list of genes to download variants, how many genes '
             'should be sent to biomart at once? Defaults to 25.'
    )
    parser.add_argument(
        '--show',
        choices=to_show,
        metavar='',
        help='For debugging only. Show chosen biomart-related data.'
             'One of following: ' + ', '.join(to_show)
    )
    parser.add_argument(
        '--profile',
        action='store_true',
        help='For debugging only.'
    )
    append_commands(parser)

    return parser


def parse_args(system_arguments):
    """Create parse, subcommands and parse given arguments list.

    Arguments list should be given in format of sys.argv."""
    subcommands = ['show', 'profile']
    parser = create_arg_parser()

    raw_args_without_filename = system_arguments[1:]

    raw_args = [
        '--' + a
        if a in subcommands
        else a
        for a in raw_args_without_filename
    ]

    return parser.parse_args(raw_args)


def show_some_variant(biomart, dataset):
    variant = None

    genes_from_patacsdb = gene_names_from_patacsdb_csv()

    for gene in genes_from_patacsdb:

        variants_by_gene = download_variants(
            biomart,
            dataset,
            [gene]
        )

        if variants_by_gene:

            gene, variants = variants_by_gene.popitem()
            variant = variants[0]
            break

    if variant:
        print(variant)
    else:
        print('No variants found in given dataset.')


@cacheable
def variants_by_gene_parsed(raw_variants_by_gene):
    from parse_variants import parse_variants_by_gene

    return parse_variants_by_gene(raw_variants_by_gene)


@cacheable
def get_all_used_transcript_ids(variants_by_gene):
    """Return all transcript identifiers that occur in the passed variants.

    variants_by_gene is expected to by a dict of variant lists grouped by genes.
    """
    transcripts_ids = set()

    for variants in variants_by_gene.itervalues():
        transcripts_ids.update(
            {
                variant.ensembl_transcript_stable_id
                for variant in variants
                }
        )

    return transcripts_ids


@cacheable
def create_cds_db(transcripts_to_load):
    return SequenceDB(
        version=ENSEMBL_VERSION,
        assembly=GRCH_VERSION,
        index_by='transcript',
        sequence_type='cds',
        restrict_to=transcripts_to_load
    )


@cacheable
def create_cdna_db(transcripts_to_load):
    return SequenceDB(
        version=ENSEMBL_VERSION,
        assembly=GRCH_VERSION,
        index_by='transcript',
        sequence_type='cdna',
        restrict_to=transcripts_to_load
    )


@cacheable
def create_dna_db():
    chromosomes = map(str, range(1, 23)) + ['X', 'Y', 'MT']
    dna_db = {}

    for chromosome in chromosomes:
        dna_db[chromosome] = FastSequenceDB(
            version=ENSEMBL_VERSION,
            assembly=GRCH_VERSION,
            sequence_type='dna',
            id_type='chromosome.' + chromosome
        )
    return dna_db


def main(args):

    snp_dataset = BiomartDataset(args.biomart, name=args.dataset)
    biomart_server = BiomartServer(args.biomart)

    show_functions = {
        'databases': biomart_server.show_databases,
        'datasets': biomart_server.show_datasets,
        'filters': snp_dataset.show_filters,
        'attributes': snp_dataset.show_attributes,
        'attributes_by_page': snp_dataset.show_attributes_by_page,
        'some_variant': lambda: show_some_variant(args.biomart, snp_dataset)
    }

    if args.show:
        func = show_functions[args.show]
        func()

    execute_commands(args)

    raw_variants_by_gene = None

    # 1. Download
    if args.download_variants:

        download_method = args.download_variants

        if download_method == 'biomart':
            genes_list = args.genes_list
            raw_variants_by_gene = download_variants.save(
                args.biomart,
                snp_dataset,
                genes_list,
                step_size=args.step_size,
                filters=args.filters
            )
        else:
            from analyses import VARIANTS_GETTERS
            raw_variants_by_gene = VARIANTS_GETTERS[download_method]()

        # transcripts have changed, some databases need reload
        transcripts_to_load = get_all_used_transcript_ids.save(raw_variants_by_gene)
        create_cds_db.save(transcripts_to_load)
        create_cdna_db.save(transcripts_to_load)

        print('Raw variants data downloaded, databases reloaded.')

    # 2. Parse
    if args.parse_variants:
        # has to be merged into one step
        download_method = args.download_variants
        if not raw_variants_by_gene:
            if download_method == 'biomart':
                raw_variants_by_gene = download_variants.load()
            else:
                from analyses import VARIANTS_GETTERS
                raw_variants_by_gene = VARIANTS_GETTERS[download_method]()

        variants_by_gene_parsed.save(raw_variants_by_gene)

        print('Variants data parsed and read to use.')

    # 3. Perform chosen analyses and generate reports
    if args.report:
        perform_analyses(args)
    else:
        print('No analyses specified.')


if __name__ == '__main__':

    parsed_args = parse_args(sys.argv)
    o.force = parsed_args.verbose
    VERBOSITY_LEVEL = parsed_args.verbose

    if parsed_args.profile:
        import profile
        profile.run('main(parsed_args)')
    else:
        main(parsed_args)
