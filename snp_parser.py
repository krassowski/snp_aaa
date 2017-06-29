#!/usr/bin/env python -u
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys
import traceback

import time

import datetime

from tqdm import tqdm
from pyfaidx import Fasta

from cache import cacheable
from commands import execute_commands, execute_subparser_commands
from commands import append_commands
from commands import append_subparsers

VERBOSITY_LEVEL = 0

GRCH_VERSION = 'GRCh37'
GRCH_SUBVERSION = '13'
ENSEMBL_VERSION = '88'
COSMIC_VERSION = '81'
DBSNP_VERSION = '150'
SPIDEX_LOCATION = 'spidex_public_noncommercial_v1.0/spidex_public_noncommercial_v1_0.tab.gz'
TRANSCRIPT_DB_PATH = 'ensembl/v' + ENSEMBL_VERSION + '/Homo_sapiens.' + GRCH_VERSION + '.cds.all.fa'
#TRANSCRIPT_DB_PATH = '/media/ramdisk/Homo_sapiens.GRCh37.cds.all.fa'


transcripts_db = Fasta(TRANSCRIPT_DB_PATH, key_function=lambda x: x.split('.')[0])


vcf_mutation_sources = {
    'COSMIC': {
        'is_alias': False,
        'path': 'cosmic/v' + COSMIC_VERSION + '/CosmicCodingMuts.vcf.gz.bgz',
        #'path': '/media/ramdisk/CosmicCodingMuts.vcf.gz.bgz',
        'given_as_positive_strand_only': True
    },
    'dbSNP': {
        'is_alias': False,
        'path': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
        GRCH_SUBVERSION + '/00-All.vcf.gz',
        'given_as_positive_strand_only': True
    },
    'ensembl': {
        'is_alias': False,
        'path': 'ensembl/v' + ENSEMBL_VERSION + '/Homo_sapiens.vcf.gz',
        #'path': '/media/ramdisk/Homo_sapiens.vcf.gz',
        'given_as_positive_strand_only': True
    },
    'ClinVar': {
        'is_alias': True,
        'aliased_vcf': 'dbSNP'
    },
    'ESP': {
        'is_alias': True,
        'aliased_vcf': 'ensembl'
    },
    'HGMD-PUBLIC': {
        'is_alias': True,
        'aliased_vcf': 'ensembl'
    },
    'PhenCode': {
        'is_alias': True,
        'aliased_vcf': 'ensembl'
    },
}


try:
    from numba import jit
except ImportError:
    print('Install numba to speed up execution')
    jit = lambda x: x


def select_poly_a_related_variants(variants):
    """Return list of variants occurring in a poly(A) track and variants which will create poly(A) track."""

    return [
        variant
        for variant in variants
        if any([
            data.has or data.will_have
            for affected_transcript in variant.affected_transcripts
            for data in affected_transcript.poly_aaa.itervalues()
        ])
    ]


def all_poly_a_variants(variants):
    poly_a_related_variants = select_poly_a_related_variants(variants.itervalues())

    for variant in poly_a_related_variants:
        yield variant


def perform_analyses(args, variants=None):
    """The main workflow happens here."""

    if not args.no_variants:
        variants = parsed_variants.load()

    from analyses import REPORTERS

    for reporter_name in args.report:
        current_reporter = REPORTERS[reporter_name]

        try:
            current_reporter(variants)

        except Exception:
            print('Analysis %s failed' % reporter_name)
            traceback.print_exc()


def create_arg_parser():
    from commands import ArgumentParserPlus
    from analyses import REPORTERS
    from variant_sources import VARIANTS_GETTERS

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
        help='Analyses to be performed; one or more from: ' + ', '.join(REPORTERS.iterkeys()),
        metavar=''
    )
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
def parsed_variants(variants):
    return variants


def main(args):

    execute_commands(args)
    execute_subparser_commands(args)

    # TODO: make it into an argument or remove later
    do_not_dump = False
    variants = None

    # 1. Download and parse
    if args.variants:

        start = time.time()

        method = args.variants

        from variant_sources import VARIANTS_GETTERS
        variants = VARIANTS_GETTERS[method](args)

        #print('Raw variants data downloaded, databases reloaded.')
        #if do_not_dump:
        #    from parse_variants import parse_variants_by_gene
        #    variants = parse_variants_by_gene(raw_variants_by_gene)
        #else:
        #    variants_by_gene_parsed.save(raw_variants_by_gene)
        parsed_variants.save(variants)

        end = time.time()
        print(
            'Variants handling finished after %s'
            %
            datetime.timedelta(seconds=end-start)
        )
        print('Variants data parsed and ready to use.')

    # 3. Perform chosen analyses and generate reports
    if args.report:
        perform_analyses(args, variants)
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

    def say(text):
        # libttspico-utils
        import os
        os.system('pico2wave -w=/tmp/x.wav "%s"' % text)
        os.system('aplay /tmp/x.wav')
        os.system('rm /tmp/x.wav')

    say('Computations have just finished!')
