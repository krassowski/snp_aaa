#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import sys
import traceback

biomart_fork_path = os.path.realpath(os.path.join(os.curdir, 'biomart'))
sys.path.insert(0, biomart_fork_path)


from biomart import BiomartServer
from tqdm import tqdm

from biomart_data import BiomartDataset
from cache import cacheable
from fasta_sequence_db import SequenceDB, FastSequenceDB
from commands import execute_commands, execute_subparser_commands
from commands import append_commands
from commands import append_subparsers

VERBOSITY_LEVEL = 0

GRCH_VERSION = 'GRCh37'
GRCH_SUBVERSION = '13'
ENSEMBL_VERSION = '75'
COSMIC_VERSION = '80'
DBSNP_VERSION = '149'
SPIDEX_LOCATION = 'spidex_public_noncommercial_v1.0/spidex_public_noncommercial_v1_0.tab.gz'


vcf_locations = {
    'COSMIC': 'cosmic/v' + COSMIC_VERSION + '/CosmicCodingMuts.vcf.gz.bgz',
    'dbSNP': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
    GRCH_SUBVERSION + '/00-All.vcf.gz',
    'ClinVar': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
    GRCH_SUBVERSION + '/00-All.vcf.gz',
    'ensembl': 'ensembl/v' + ENSEMBL_VERSION + '/Homo_sapiens.vcf.gz'
}


try:
    from numba import jit
except ImportError:
    print('Install numba to speed up execution')
    jit = lambda x: x


def select_poly_a_related_variants(variants):
    """Return list oof variants occurring in a poly(A) track and variants which will create poly(A) track."""
    from parse_variants import analyze_poly_a

    return [
        variant
        for variant in map(analyze_poly_a, variants)
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

    if not args.no_variants:
        variants_by_gene = variants_by_gene_parsed.load()
    else:
        variants_by_gene = None

    from analyses import REPORTERS

    for reporter_name in args.report:
        current_reporter = REPORTERS[reporter_name]

        try:
            current_reporter(variants_by_gene)

        except Exception:
            traceback.print_exc()


def create_arg_parser():
    from commands import ArgumentParserPlus
    from analyses import REPORTERS
    from variant_sources import VARIANTS_GETTERS

    to_show = [
        'databases', 'datasets', 'filters', 'attributes',
        'attributes_by_page',
    ]

    parser = ArgumentParserPlus(description=(
        'Retrieve all data about variants from genes belonging to a predefined '
        'or user-chosen list of genes (see --genes_list) [including list of all '
        'human genes] and perform one or multiple of available analyses '
        '(see --report).\n'
        'Once variants are downloaded and pre-parsed, they will be stored '
        'in cache until another list of variants will be specified with '
        'options: --variants.'
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
        '--verbose',
        action='store_true',
        help='increase output verbosity'
    )
    parser.add_argument(
        '-n',
        '--no_variants',
        action='store_true',
        help='Use to run standalone analysis without feeding it with variants'
    )
    parser.add_argument(
        '-v', '--variants',
        choices=VARIANTS_GETTERS.keys(),
        default=None,
        help='How should the variants be retrieved? Choices are: ' +
             ', '.join([
                '{getter_name} ({description})'.format(
                    getter_name=name,
                    description=func.__doc__
                )
                for name, func in VARIANTS_GETTERS.items()
            ])
    )
    parser.add_argument(
        '--show',
        choices=to_show,
        metavar='',
        help='For debugging only. Show chosen biomart-related data.'
             ' One of following: ' + ', '.join(to_show)
    )
    parser.add_argument(
        '--profile',
        action='store_true',
        help='For debugging only.'
    )
    append_commands(parser)
    append_subparsers(parser)

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
                transcript.ensembl_id
                for variant in variants
                for transcript in variant.affected_transcripts
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

    if args.variants == 'biomart':
        snp_dataset = BiomartDataset(args.biomart, name=args.dataset)
        biomart_server = BiomartServer(args.biomart)

        show_functions = {
            'databases': biomart_server.show_databases,
            'datasets': biomart_server.show_datasets,
            'filters': snp_dataset.show_filters,
            'attributes': snp_dataset.show_attributes,
            'attributes_by_page': snp_dataset.show_attributes_by_page,
        }

        if args.show:
            func = show_functions[args.show]
            func()

    execute_commands(args)
    execute_subparser_commands(args)

    # 1. Download and parse
    if args.variants:

        method = args.variants

        from variant_sources import VARIANTS_GETTERS
        raw_variants_by_gene = VARIANTS_GETTERS[method](args)

        # transcripts have changed, some databases need reload
        transcripts_to_load = get_all_used_transcript_ids.save(raw_variants_by_gene)
        create_cds_db.save(transcripts_to_load)
        create_cdna_db.save(transcripts_to_load)

        print('Raw variants data downloaded, databases reloaded.')

        variants_by_gene_parsed.save(raw_variants_by_gene)

        print('Variants data parsed and ready to use.')

    # 3. Perform chosen analyses and generate reports
    if args.report:
        perform_analyses(args)
    else:
        print('No analyses specified.')


if __name__ == '__main__':

    parsed_args = parse_args(sys.argv)
    VERBOSITY_LEVEL = parsed_args.verbose

    if parsed_args.profile:
        import profile
        profile.run('main(parsed_args)')
    else:
        main(parsed_args)
