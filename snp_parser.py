#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import traceback
from argparse import ArgumentParser

from multiprocessing import freeze_support

from analyses import analyses
from cache import cacheable
from commands import append_commands
from commands import append_subparsers
from commands import execute_commands, execute_subparser_commands
from helpers import execution_time
from variant_sources import sources


@cacheable
def variants_cache(variants):
    return variants


def perform_analyses(args, variants=None):
    """The main workflow happens here."""

    if not args.no_variants:
        variants = variants_cache.load()

    for reporter_name in args.report:
        current_reporter = analyses[reporter_name]

        try:
            current_reporter(variants)
        except Exception:
            print('Analysis %s failed' % reporter_name)
            traceback.print_exc()


def create_arg_parser():

    parser = ArgumentParser(description=(
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
        choices=analyses.keys(),
        default=[],
        help='Analyses to be performed; one or more from: ' + ', '.join(analyses.keys()),
        metavar=''
    )
    parser.add_argument(
        '-n',
        '--no_variants',
        action='store_true',
        help='Use to run standalone analysis without feeding it with variants'
    )
    parser.add_argument(
        '-v', '--variants',
        choices=sources.keys(),
        default=None,
        help='How variants should be retrieved? Choices are: ' +
             ', '.join([
                '{source_name} ({description})'.format(
                    source_name=name,
                    description=func.__doc__
                )
                for name, func in sources.items()
                ]
             )
    )
    append_commands(parser)
    append_subparsers(parser)

    return parser


def main(args):

    execute_commands(args)
    execute_subparser_commands(args)

    variants = None

    # Download and parse
    if args.variants:

        with execution_time() as time:
            method = args.variants

            variants = sources[method](args)

            variants_cache.save(variants)

        print('Variants handling finished after %s' % time.elapsed)
        print('Variants data parsed and ready to use.')

    # Perform chosen analyses and generate reports
    if args.report:
        perform_analyses(args, variants)
    else:
        print('No analyses specified.')

if __name__ == '__main__':
    freeze_support()
    parser = create_arg_parser()
    parsed_args = parser.parse_args()

    main(parsed_args)

