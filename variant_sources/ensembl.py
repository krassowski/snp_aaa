from __future__ import print_function
import gzip
import json
import subprocess
from collections import defaultdict, OrderedDict
import gc
from contextlib import contextmanager
from tqdm import tqdm
from cache import args_aware_cacheable
from commands import SourceSubparser
from variant import Variant, AffectedTranscript
from variant_sources import variants_getter
from variant_sources.biomart import gene_names_from_patacsdb_csv


@contextmanager
def fast_gzip_read(file_name, single_thread=False):
    command = 'zcat %s' if single_thread else 'unpigz -p 8 -c %s'
    print(command % file_name)
    p = subprocess.Popen(
        command % file_name,
        shell=True,
        bufsize=100000000,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    yield p.stdout


@args_aware_cacheable
def count_lines(file_name, single_thread=False):
    command = 'zcat %s | wc -l' if single_thread else 'unpigz -p 8 -c %s | wc -l'

    out = subprocess.Popen(
        command % file_name,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
     ).communicate()[0]

    return int(out.partition(b' ')[0])


ensembl_args = SourceSubparser(
    'ensembl',
    help='Arguments for raw Ensembl variants loading'
)

ensembl_args.add_command(
    '--early_selection',
    choices=['spidex', 'full_aaa'],
    default='spidex',
    help=(
        'Loading all variants takes ages. "early_selection" option specifies analysis for which the variants will '
        'be used so loader can choose which variants should be rejected early. '
        'Default: spidex. Choices: spidex, full_aaa.'
    )
)

"""
ensembl_args.add_command(
    '--filters',
    type=json.loads,
    default={},
    help='Additional variants filters to be used when querying biomart'
)

ensembl_args.add_command(
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
"""


@variants_getter
def ensembl(args):
    """Load variants from Ensembl raw database files.
    Only variants belonging to coding sequences will be loaded.

    """

    # specified with --genes_list and filtered with --filters.
    # gene_names = args.genes_list
    # filters = args.filters

    types = OrderedDict(
        (
            ('transcript_variation_id', int),
            ('variation_feature_id', int),
            ('feature_stable_id', str),
            ('allele_string', str),
            ('somatic', bool),
            ('consequence_types', str),
            ('cds_start', int),
            ('cds_end', int),
            ('cdna_start', int),
            ('cdna_end', int),
            ('translation_start', int),
            ('translation_end', int),
            ('distance_to_transcript', int),
            ('codon_allele_string', str),
            ('pep_allele_string', str),
            ('hgvs_genomic', str),
            ('hgvs_transcript', str),
            ('hgvs_protein', str),
            ('polyphen_prediction', str),
            ('polyphen_score', float),
            ('sift_prediction', str),
            ('sift_score', float),
            ('display', bool)
        )
    )
    loc = 'ensembl/v88/GRCh37/variation_database/'

    by_id = {}
    keys = types.keys()
    somatic_pos = keys.index('somatic')
    consequence_pos = keys.index('consequence_types')

    from snp_parser import jit

    @jit
    def int_or_none(x):
        if x == '\\N':
            return None
        return int(x)

    transcript_strand = {}
    filename = loc + 'transcript.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            transcript_strand[data[14]] = int(data[6])

    accepted_consequences = {
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
        'stop_lost',
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
    }

    @jit
    def good_consequence(data):
        return not accepted_consequences.isdisjoint(data[consequence_pos].split(','))

    @jit
    def is_spidex_poly_aaa_viable(alleles):
        return (
            len(alleles[0]) == 1 and
            (
                (
                    len(alleles[1]) == 1 and
                    'A' in alleles
                ) or
                alleles[1] == 'COSMIC_MUTATION'
            )
        )

    @jit
    def is_poly_aaa_viable(alleles):
        has_a = False
        for allele in alleles:
            if 'A' in allele:
                has_a = True
                break

        return (
            has_a
            # or alleles[1] == 'COSMIC_MUTATION'    # this will be included in has_a ;)
        )

    if args.early_selection == 'spidex':
        is_viable = is_spidex_poly_aaa_viable
    else:
        is_viable = is_poly_aaa_viable

    # TODO: LRG count them, ignore mismatches or ignore them at all
    # TODO vcf parsers: ise ensembl needed?

    filename = loc + 'transcript_variation.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            #if data[somatic_pos] == '1' and good_consequence(data):
            if good_consequence(data):
                alleles = data[3].split('/')

                if not is_viable(alleles):
                    continue

                ref = alleles[0]
                variant_id = int(data[1])

                transcript = AffectedTranscript(
                    cds_start=int_or_none(data[6]),
                    cds_end=int_or_none(data[7]),
                    # cdna_start=int_or_none(data[8]),
                    # cdna_end=int_or_none(data[9]),
                    strand=transcript_strand[data[2]],
                    ensembl_id=data[2]
                )

                if variant_id in by_id:
                    v = by_id[variant_id]

                    if v.allele_1 != ref:
                        print(
                            'Reference does not match: %s vs %s'
                            'for variant %s, with transcripts: %s %s'
                            %
                            (
                                v.allele_1,
                                ref,
                                variant_id,
                                data[2],
                                list(v.affected_transcripts)[0].ensembl_id
                            )
                        )
                        continue

                    v.affected_transcripts.add(transcript)
                    # if v.ensembl_transcript_stable_id == data[2]:
                    #    v
                else:
                    v = Variant(
                        allele_1=ref,
                        affected_transcripts={transcript}
                    )
                    by_id[variant_id] = v
    gc.collect()

    print('Accepted:', len(by_id))

    seq_region = {}
    # seq_region_id, name, cord_system_fk
    filename = loc + 'seq_region.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            seq_region[int(data[0])] = data[1]

    sources = {}
    # 'source_id', 'name', 'version', 'description', 'url', 'type', 'somatic_status', 'data_types'
    filename = loc + 'source.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            sources[int(data[0])] = data[1]
    gc.collect()

    feature_table = [
        'variation_feature_id',
        'seq_region_id',
        'seq_region_start',
        'seq_region_end',
        'seq_region_strand',
        'variation_id',
        'allele_string',
        'variation_name',
        'map_weight',
        'flags',
        'source_id',
        'consequence_types',
        'variation_set_id',
        'class_attrib_id',
        'somatic',
        'minor_allele',
        'minor_allele_freq',
        'minor_allele_count',
        'alignment_quality',
        'evidence_attribs',
        'clinical_significance',
        'display',
    ]

    filename = loc + 'variation_feature.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            if int(data[0]) in by_id:
                v = by_id[int(data[0])]
                v.chr_name = seq_region[int(data[1])]
                v.chrom_start = int(data[2])
                v.chrom_end = int(data[3])
                v.chrom_strand = data[4]
                v.refsnp_id = data[7]
                v.refsnp_source = sources[int(data[10])]
                # if v.refsnp_source != 'COSMIC':
                #    print(v)

    gc.collect()
    by_name = defaultdict(list)
    no_name = 0
    no_u = 0

    for v in tqdm(by_id.itervalues(), total=len(by_id)):
        if not hasattr(v, 'refsnp_id'):
            no_name += 1
            continue
        if v.refsnp_id in by_name:
            by_name[v.refsnp_id].affected_transcripts.update(v.affected_transcripts)
            no_u += 1
            #print('Compare:')
            #print(by_name[v.refsnp_id])
            #print(v)
        else:
            by_name[v.refsnp_id] = v
            # transcript = v.extract_transcript()
    print('Noname: %s' % no_name)
    print('Nouniqe names: %s' % no_u)
    print('Remained grouped by name: %s' % len(by_name))

    #del by_id
    gc.collect()
    #print('Press enter to continue')
    #x = raw_input()

    # organise variants in groups (not by genes but randomly) for better use of multi-threading
    out = {}
    c = 0
    group = []
    ids = []
    for v_id, v in by_name.iteritems():
        c += 1
        if c % 5000 == 0:
            out[','.join(ids)] = group
            group = []
            ids = []
        group.append(v)
        ids.append(v_id)

    out[','.join(ids)] = group
    # out = {v_id: [v] for v_id, v in by_name.iteritems()}
    print('Enesmbl variants ready for parsing')
    return out
    # THERE ARE MANY GENES FOR A SINGLE VARIANT, HENCE I DO NOT ORGANISE VARIANTS BY GENES
