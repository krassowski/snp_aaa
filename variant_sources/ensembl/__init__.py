from __future__ import print_function

import gc
from collections import defaultdict

from tqdm import tqdm

from commands import SourceSubparser
from multiprocess import fast_gzip_read, count_lines, get_manager
from variant_sources import variants_source
from variant_sources.ensembl.load_transcript_variant_pairs import load_poly_a_transcript_variant_pairs
from variant_sources.ensembl.load_varaints_details import load_variants_details
from variant_sources.ensembl.poly_aaa import get_poly_aaa


def transcript_names_from_patacsdb_csv(filename='patacsdb_all.csv'):
    """
    Load stable ensembl transcript ids from given csv file.
    The default file was exported from homo sapiens dataset
    of patacsdb database (sysbio.ibb.waw.pl/patacsdb).

    The identifiers are expected to be located in third column of each row.
    """
    transcripts = set()
    with open(filename) as f:
        header = next(f)     # skip header
        assert header == '"Protein name","Gene id","Transcript id","Location (%)",Sequence'
        for line in f:
            data = line.strip().split(',')
            transcripts.add(data[2])

    return list(transcripts)

ensembl_args = SourceSubparser(
    'ensembl',
    help='Arguments for raw Ensembl variants loading'
)

ensembl_args.add_command(
    '--location',
    default='ensembl/v88/GRCh37/variation_database/'
)


def set_maker(string):
    return set(string.split(','))


ensembl_args.add_command(
    '--consequences',
    default={
        'coding_sequence_variant',
        'synonymous_variant',
        'stop_retained_variant',
        'protein_altering_variant',
        'inframe_variant',
        'incomplete_terminal_codon_variant',
        'missense_variant',
        'non_conservative_missense_variant',
        'conservative_missense_variant'
        'stop_gained',
        'initiator_codon_variant',
        'inframe_indel',
        'inframe_insertion',
        'conservative_inframe_insertion',
        'disruptive_inframe_insertion',
        'inframe_deletion',
        'disruptive_inframe_deletion',
        'conservative_inframe_deletion',
        'frameshift_variant',
        'plus_1_frameshift_variant',
        'minus_1_frameshift_variant',
        'frameshift_elongation',
        'frame_restoring_variant',
        'plus_2_frameshift_variant',
        'frameshift_truncation',
        'minus_2_frameshift_variant',
        'terminator_codon_variant',
        'incomplete_terminal_codon_variant',
        'stop_lost',
        'stop_retained_variant',
    },
    type=set_maker
)


@ensembl_args.command('--transcript_stats', action='store_true')
def transcript_stats(value, args):
    if not value:
        return

    print('Generating transcripts statistics')

    loc = args.location
    filename = loc + 'transcript_variation.txt.gz'
    all = float(count_lines(filename))
    enst = count_lines(filename, grep='ENST')
    lrg = count_lines(filename, grep='LRG_')
    print('All transcripts: %s ' % all)
    print('ENST transcripts: %s (%s percent)' % (enst, enst/all))
    print('LRG transcripts: %s (%s percent)' % (lrg, lrg/all))
    print('Other transcripts: %s (%s percent)' % (all-enst-lrg, (all-lrg-enst)/all))


ensembl_args.add_command(
    '--transcripts',
    nargs='+',
    type=str,
    default='all',
    help=(
        'list of human genes from which variants should be extracted for use in '
        'available analyses. By default all human transcripts are used. '
        'To restrict to PATACSDB transcript use "patacsdb".'
    )
)


ensembl_args.add_command(
    '--not_only_poly_a',
    action='store_true',
    help=(
        'include all variants rather than only those which are poly_a related'
    )
)


def load_variation_sources(loc):
    sources = {}
    # 'source_id', 'name', 'version', 'description', 'url', 'type', 'somatic_status', 'data_types'
    filename = loc + 'source.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            sources[int(data[0])] = data[1]
    gc.collect()
    return sources


def load_chromosome_and_region_names(loc):
    seq_region = {}
    # seq_region_id, name, cord_system_fk
    filename = loc + 'seq_region.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            seq_region[int(data[0])] = data[1]
    return seq_region


def load_transcript_strands(loc):
    transcript_strand = {}
    filename = loc + 'transcript.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            transcript_strand[data[14]] = int(data[6])
    return transcript_strand


def aggregate_by_variant(transcript_variant_pairs):
    print('Aggregating variant transcript pairs')
    variants_by_id = defaultdict(list)
    for internal_variant_id, transcript, alleles in tqdm(transcript_variant_pairs):
        variants_by_id[internal_variant_id].append((transcript, alleles))

    return variants_by_id


@variants_source
def ensembl(args):
    """Load variants from Ensembl raw database files.
    Only variants belonging to coding sequences will be loaded.
    """
    print(args)

    path = args.location
    only_poly_a = not args.not_only_poly_a

    if args.transcripts == 'all':
        only_transcripts = None
    elif args.transcripts == 'patacsdb':
        only_transcripts = transcript_names_from_patacsdb_csv()
    else:
        only_transcripts = args.transcripts

    transcript_strand = load_transcript_strands(path)
    seq_region = load_chromosome_and_region_names(path)
    sources = load_variation_sources(path)

    # multiple threads:
    transcript_variant_pairs = load_poly_a_transcript_variant_pairs(
        path, transcript_strand, args.consequences, only_transcripts, only_poly_a
    )
    accepted_transcript_variant_pairs, to_check_transcript_variant_pairs = transcript_variant_pairs

    # those two can run in two separate processes but both have to be single-threaded.
    accepted_transcripts_by_variant_id = aggregate_by_variant(accepted_transcript_variant_pairs)
    to_check_transcripts_by_variant_id = aggregate_by_variant(to_check_transcript_variant_pairs)

    # having a dict passed to each instance can use a lot of memory; having a proxy slows down access
    # using proxy dict for the bigger dict and normal dict for the smaller one seems like a good trade-of
    # there are usually more mutations to check than those which are accepted:
    manager = get_manager()
    to_check_transcripts_by_variant_id = manager.dict(to_check_transcripts_by_variant_id)

    del accepted_transcript_variant_pairs
    del to_check_transcript_variant_pairs
    del transcript_variant_pairs

    # multiple threads:
    variants_by_id = load_variants_details(
        path,
        seq_region, sources,
        accepted_transcripts_by_variant_id, to_check_transcripts_by_variant_id,
        only_poly_a
    )

    print('Accepted:', len(variants_by_id))

    variants_by_name = group_variants_by_snp_id(variants_by_id)

    print('Enesmbl variants ready for parsing')
    return variants_by_name


def group_variants_by_snp_id(by_id):
    by_name = defaultdict(list)
    no_name = 0
    no_u = 0
    for raw_id in tqdm(by_id.keys(), total=len(by_id)):
        v = by_id[raw_id]
        if not v.snp_id:
            no_name += 1
            continue
        if v.snp_id in by_name:
            by_name[v.snp_id].affected_transcripts.extend(v.affected_transcripts)
            no_u += 1
        else:
            by_name[v.snp_id] = v
    if no_name:
        print('Variant-transcript pairs without SNP id: %s' % no_name)
        print('Are you analysing head?')
    if no_u:
        print('SNP ids having more than one variant-transcript pair: %s' % no_u)

    print('Number of unique SNP_id-based variant groups: %s' % len(by_name))
    gc.collect()
    return by_name


