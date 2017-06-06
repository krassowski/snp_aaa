from __future__ import print_function
import subprocess
from collections import defaultdict, OrderedDict
import gc
from contextlib import contextmanager
from multiprocessing import Manager, Process

import itertools
from tqdm import tqdm
from cache import args_aware_cacheable
from commands import SourceSubparser
from parse_variants import OFFSET, complement
from snp_parser import transcripts, vcf_mutation_sources
from variant import Variant, AffectedTranscript, PolyAAAData
from variant_sources import variants_getter
from snp_parser import jit
from vcf_parser import VariantCallFormatParser, ParsingError
from poly_a import poly_a


def show_context(seq, start=OFFSET, end=-OFFSET):
    print(seq[:start] + '>' + seq[start:end] + '<' + seq[end:])


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


@contextmanager
def fast_gzip_read(file_name, single_thread=False):
    command = 'zcat %s' if single_thread else 'unpigz -p 4 -c %s'
    print(command % file_name)
    p = subprocess.Popen(
        (command % file_name).split(' '),
        #shell=True,
        #bufsize=500000,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    yield p.stdout


@args_aware_cacheable
def count_lines(file_name, single_thread=False, grep=None):

    command = 'zcat %s' if single_thread else 'unpigz -p 6 -c %s'
    if grep:
        command += '| grep %s' % grep
    command += ' | wc -l'

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
    choices=['spidex_poly_aaa', 'all_potential_poly_aaa'],
    default='spidex_poly_aaa',
    help=(
        'Loading all variants takes ages. "early_selection" option specifies analysis for which the variants will '
        'be used so loader can choose which variants should be rejected early. '
        'Default: spidex. Choices: spidex, all_potential_poly_aaa.'
    )
)

ensembl_args.add_command(
    '--location',
    #default='/media/ramdisk/head/'
    default='ensembl/v88/GRCh37/variation_database/head/'
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


"""
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


@jit
def get_sequence(transcript, raw_start, raw_end):

    if not (raw_start and raw_end):
        # print('No sequence coordinates for %s transcript' % transcript.ensembl_id)
        return

    if transcript not in transcripts:
        return

    if raw_end - raw_start > 2 * OFFSET:
        # print('Skipping transcript %s of variant %s: too wide variation' % (transcript.ensembl_id, variant.snp_id))
        return

    whole_seq = str(transcripts[transcript])

    if whole_seq:
        start = raw_start - 1
        end = raw_end

        cut_from = start - OFFSET
        cut_to = end + OFFSET

        seq = whole_seq[max(cut_from, 0):max(cut_to, 0)]

        # if we are on the edge of sequence
        if cut_from < 0:
            seq = '-' * (-cut_from) + seq
        if cut_to > len(whole_seq):
            seq += '-' * (cut_to - len(whole_seq))

        assert len(seq) == OFFSET * 2 + raw_end - raw_start + 1
    else:
        return

    return seq


@jit
def get_poly_a(sequence, alternative_alleles):

    has_aaa, before_len = poly_a(
        sequence,
        OFFSET,
        len(sequence) - OFFSET
    )

    poly_aaa = defaultdict(PolyAAAData)

    for alt in alternative_alleles:

        mutated_seq = sequence[:OFFSET] + alt + sequence[-OFFSET:]

        will_have, after_len = poly_a(
            mutated_seq,
            OFFSET,
            len(mutated_seq) - OFFSET
        )

        poly_aaa[alt].has = has_aaa
        poly_aaa[alt].will_have = will_have
        poly_aaa[alt].before = before_len
        poly_aaa[alt].after = after_len

    return poly_aaa


class AlleleMismatch(Exception):
    pass


class ToFewAlleles(Exception):
    pass


def parse_alleles(alleles, transcript, convert_all=True):
    sequence_ref_allele = transcript.sequence[OFFSET:-OFFSET]

    if not alleles:
        print('No alleles for:')
        return

    #if len(alleles) == 1:
    #    raise ToFewAlleles()

    if len(alleles) == 2:
        allele = alleles[1]

        if allele.startswith('(') and allele.endswith(')'):
            alleles = parse_short_tandem_repeat(allele, alleles)

    alleles = ['' if a == '-' else a for a in alleles]

    if transcript.strand == -1:
        if convert_all:
            alleles = [
                complement(alt)[::-1] for alt in alleles
            ]
        else:
            alleles = [alleles[0]] + [
                complement(alt)[::-1] for alt in alleles[1:]
            ]

    if sequence_ref_allele != alleles[0]:
        raise AlleleMismatch(sequence_ref_allele, alleles[0])

    return alleles


def parse_short_tandem_repeat(allele, alleles):
    """
    Reference: https://www.ncbi.nlm.nih.gov/variation/hgvs/

    So we have such a case:
    59521612	27517	138664867	138664881	-1	60499308	(AGCTGCGGCTGCAGC(3))	rs387906322	1	\N	1	coding_sequence_variant,upstream_gene_variant,regulatory_region_variant	1,23,24,30,31,32,41	18	0	\N	\N	\N	1	418	pathogenic	1
    where the variant is a short tandem repeat.
    here is corresponding page: https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=rs387906322

    For repeats where there is a range given, an average repeat count (round down) will be used:
    50313977	27522	4680020	4680043	1	51297695	(CAGGGCGGTGGTGGCTGGGGGCAG(6_13))	rs367543047	1	\N	1	coding_sequence_variant	1,23,30,31,32,41	18	0	\N\N	\N	1	371,418	pathogenic	1
    between 6 and 13 copies http://varnomen.hgvs.org/recommendations/DNA/variant/repeated/
    so (6 + 13) / 2 ~= 10

    Args:
        allele:
        alleles:

    Returns:
    """
    allele = allele[1:-1]
    if '(' in allele:
        ref, repeat_count = allele.split('(')
        repeat_count = repeat_count[:-1]
        if '_' in repeat_count:
            min_count, max_count = map(int, repeat_count.split('_'))
            repeat_count = (max_count + min_count) / 2
        else:
            repeat_count = int(repeat_count)
        alt = ref * repeat_count
        assert alleles[0] == ref
        alleles = [ref, alt]
        return alleles
    elif ' BP DELETED' in allele:
        cnt = len(alleles[0])
        if ('%s BP DELETED' % cnt) == allele:
            return [allele[0], '']
        else:
            raise ParsingError('Not a STR, neither known exception', alleles, allele)
    else:
        raise ParsingError('Failed to parse STR', alleles, allele)


# Warning: adding jit causes wrong behaviour!
#@jit
def analyze_poly_aaa(transcript, alleles):
    """Add poly_aaa to AffectedTranscript if possible.
    Returns True if poly_aaa has been detected, False-considered value otherwise."""
    transcript.poly_aaa = get_poly_a(transcript.sequence, alleles[1:])

    for poly_a in transcript.poly_aaa.values():
        if poly_a.has or poly_a.will_have:
            #print('True')
            return True


vcf_parser = VariantCallFormatParser(vcf_mutation_sources)


def get_alt_alleles(variant, strand, convert_to_strand):

    vcf_data = vcf_parser.get_by_variant(variant)
    if len(vcf_data) > 1:
        print(
            'VCF data contains more than one record matching %s '
            'variant: %s ' % (variant.snp_id, vcf_data)
        )
        print('Only the first record will be parsed!')

    elif len(vcf_data) == 0:
        raise ParsingError('No VCF data for %s.' % variant.snp_id)

    analysed_vcf_record = vcf_data[0]
    pos, ref, alts = vcf_parser.parse(analysed_vcf_record, variant.source, strand, convert_to_strand)
    # gene = vcf_parser.get_gene(analysed_vcf_record)
    return pos, ref, alts


class IncorrectAllele(Exception):
    pass


@jit
def int_or_none(x):
    if x == '\\N':
        return None
    return int(x)


@variants_getter
def ensembl(args):
    """Load variants from Ensembl raw database files.
    Only variants belonging to coding sequences will be loaded.

    """

    print(args)

    loc = args.location
    accepted_consequences = args.consequences

    # specified with --genes_list and filtered with --filters.
    # gene_names = args.genes_list
    # filters = args.filters

    transcript_strand = {}
    filename = loc + 'transcript.txt.gz'
    with fast_gzip_read(filename) as f:
        for line in tqdm(f, total=count_lines(filename)):
            data = line.split('\t')
            transcript_strand[data[14]] = int(data[6])

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
    keys = types.keys()
    somatic_pos = keys.index('somatic')
    consequence_pos = keys.index('consequence_types')

    @jit
    def good_consequence(data):
        # to restrict to somatic mutations, use: and data[somatic_pos] == '1'
        return not accepted_consequences.isdisjoint(data[consequence_pos].split(','))

    filename = loc + 'transcript_variation.txt.gz'
    by_id, to_check_later = load_variants_and_transcripts_if_poly_aaa(
        filename, good_consequence, transcript_strand
    )
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

    filename = loc + 'variation_feature.txt.gz'
    load_variants_details(by_id, filename, seq_region, sources, to_check_later)
    gc.collect()

    by_name = group_variants_by_snp_id(by_id)
    gc.collect()

    out = group_variants_in_chunks(by_name)

    print('Enesmbl variants ready for parsing')
    return out
    # THERE ARE MANY GENES FOR A SINGLE VARIANT, HENCE I DO NOT ORGANISE VARIANTS BY GENES


def save_if_poly_a_related(by_id, transcript, parsed_alleles, internal_variant_id):
    detected = analyze_poly_aaa(transcript, parsed_alleles)

    if not detected:
        return

    print('Poly(A) detected in:', transcript.sequence, internal_variant_id)

    ref = parsed_alleles[0]
    # TODO: reintroduce this check later in sanity check
    #"""
    v = by_id.get(internal_variant_id, None)

    if v:
        if not v.ref:
            v.ref = ref
        elif v.ref != ref:
            print(
                'Reference does not match: %s vs %s '
                'for variant %s, with transcripts: %s %s'
                %
                (
                    v.ref,
                    ref,
                    internal_variant_id,
                    transcript.ensembl_id,
                    ', '.join(t.ensembl_id for t in v.affected_transcripts)
                )
            )
            return

        v.affected_transcripts.add(transcript)
        by_id[internal_variant_id] = v
    else:
        #"""
        by_id[internal_variant_id] = Variant(
            ref=ref,
            affected_transcripts={transcript}
        )
    return True


def grouper(iterable, chunk_size, fill_value=None):
    args = [iter(iterable)] * chunk_size
    return itertools.izip_longest(fillvalue=fill_value, *args)


num_workers = 6
manager = Manager()


def load_variants_and_transcripts_if_poly_aaa(
        filename, good_consequence, transcript_strand):

    def do_work(in_queue, by_id, to_check_later):
        while True:
            lines = in_queue.get()

            if lines is None:
                in_queue.task_done()
                return

            for line in lines:

                if line is None:
                    in_queue.task_done()
                    return

                data = line.split('\t')

                # ignore 'LRG' transcripts at all
                if data[2].startswith('LRG_'):
                    continue

                if not good_consequence(data):
                    continue

                alleles = data[3].split('/')

                ensembl_id = data[2]

                transcript = AffectedTranscript(
                    cds_start=int_or_none(data[6]),
                    cds_end=int_or_none(data[7]),
                    strand=transcript_strand[ensembl_id],
                    ensembl_id=ensembl_id
                )

                sequence = get_sequence(transcript.ensembl_id, transcript.cds_start, transcript.cds_end)

                if not sequence:
                    continue

                transcript.sequence = sequence

                internal_variant_id = int(data[1])

                if len(alleles) < 1:
                    # 152791410	27516	66766360	66766362	1	153751119	(GGC)18	rs869109080	1	\N	1	non_coding_transcript_variant,non_coding_transcript_exon_variant,NMD_transcript_variant,inframe_insertion,regulatory_region_variant	1	16	0	\N	\N	\N	1	\N	\N	0
                    continue

                # TODO: why some PhenCode variations ends with 'puh'?
                if alleles[1] in ('COSMIC_MUTATION', 'HGMD_MUTATION', 'PhenCode_variation') or alleles[1].endswith('puh'):
                    if internal_variant_id not in to_check_later:
                        to_check_later[internal_variant_id] = list()
                    to_check_later[internal_variant_id].append(
                        (transcript, alleles)
                    )
                    continue

                try:
                    alleles = parse_alleles(alleles, transcript)
                except AlleleMismatch:
                    # it happens that alleles are not in following order: ref/(alts)
                    # it is rare, but there is no documentations which would forbid such annotation. An example:
                    # https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=rs11543174
                    # because checking which allele is reference allele in VCF file takes quite a long time, I cannot
                    # am reluctant to check this for each variant. Instead, checking it only for rare cases like this
                    # makes whole process way much faster and still reliable.
                    if internal_variant_id not in to_check_later:
                        to_check_later[internal_variant_id] = list()
                    to_check_later[internal_variant_id].append(
                        (transcript, alleles)
                    )
                    continue
                except (KeyError, ParsingError) as e:
                    print(e)
                    print(alleles)
                    print(internal_variant_id)
                    print(transcript)
                    show_context(transcript.sequence)
                    continue

                save_if_poly_a_related(by_id, transcript, alleles, internal_variant_id)

    by_id = manager.dict()
    to_check_later = manager.dict()
    work = manager.Queue(num_workers)

    pool = []
    for i in xrange(num_workers):
        p = Process(target=do_work, args=(work, by_id, to_check_later))
        p.start()
        pool.append(p)

    chunk_size = 100000
    with fast_gzip_read(filename) as f:
        iterator = grouper(f, chunk_size)
        for i in tqdm(iterator, total=count_lines(filename) / chunk_size):
            work.put(i)

    print('Work distributed')

    for i in xrange(num_workers):
        work.put(None)

    for p in pool:
        p.join()

    print('Joined all')

    return by_id, to_check_later


def group_variants_in_chunks(by_name, chunk_size=15000):
    """Organise variants in groups (not by genes but randomly) for better use of multi-threading"""
    out = {}
    c = 0
    group = []
    ids = []
    for v_id, v in by_name.iteritems():
        c += 1
        if c % chunk_size == 0:
            out[','.join(ids)] = group
            group = []
            ids = []
        group.append(v)
        ids.append(v_id)
    out[','.join(ids)] = group
    return out


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
            by_name[v.snp_id].affected_transcripts.update(v.affected_transcripts)
            no_u += 1
        else:
            by_name[v.snp_id] = v
    if no_name:
        print('Variant-transcript pairs without SNP id: %s' % no_name)
        print('Are you analysing head?')
    if no_u:
        print('SNP ids having more than one variant-transcript pair: %s' % no_u)

    print('Number of unique SNP_id-based variant groups: %s' % len(by_name))
    return by_name


def load_variants_details(by_id, filename, seq_region, sources, to_check_later):
    """
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

    Args:
        by_id:
        filename:
        seq_region:
        sources:
        to_check_later:

    Returns:
    """

    def do_work(in_queue, by_id, to_check_later):
        while True:
            lines = in_queue.get()

            if lines is None:
                in_queue.task_done()
                return

            for line in lines:

                if line is None:
                    in_queue.task_done()
                    return

                data = line.split('\t')
                variant_id = int(data[0])

                in_by_id = variant_id in by_id
                in_to_check_later = variant_id in to_check_later

                if in_to_check_later and (not in_by_id):
                    by_id[variant_id] = Variant()
                    in_by_id = True

                if in_by_id:
                    variant = by_id[variant_id]

                    variant.chr_name = seq_region[int(data[1])]
                    variant.chr_start = int(data[2])
                    variant.chr_end = int(data[3])
                    variant.chr_strand = data[4]
                    variant.snp_id = data[7]
                    variant.source = sources[int(data[10])]

                    by_id[variant_id] = variant

                if in_to_check_later:
                    saved = False

                    checked_transcripts = set()

                    for transcript, alleles in to_check_later[variant_id]:
                        checked_transcripts.add(transcript)
                        try:
                            pos, ref, alts = get_alt_alleles(variant, transcript.strand, convert_to_strand=False)
                        except ParsingError as e:
                            if variant.source == 'COSMIC':
                                # TODO: more specific exception
                                # some COSMIC mutations are mapped twice to ensembl but only one
                                # position is given in Cosmic VCF file. In such cases, I stick with
                                # Cosmic VCF file. There is no need to print out such cases
                                continue
                            else:
                                print('ParsingError:', e.message, e.args)
                                print(variant_id, variant.snp_id, alleles)
                            continue

                        if any('.' in alt for alt in alts):
                            if variant.source == 'HGMD-PUBLIC':
                                pass
                                # TODO: count this
                                # print('Skipping variant from HGMD without complete allele data: %s' % v.snp_id)
                            else:
                                raise ParsingError('IncorrectAllele, with dot, not from HGMD')
                            continue

                        raw_alleles = [alleles[0]] + alts

                        try:
                            alleles = parse_alleles(raw_alleles, transcript)
                        except AlleleMismatch as e:
                            if variant.source == 'COSMIC':
                                # sometimes ensembl has wrongly stated ref allele too.
                                print('Variant %s from COSMIC had wrongly stated ref alleles' % variant.snp_id)
                                try:
                                    alleles = parse_alleles(raw_alleles, transcript, convert_all=False)
                                except AlleleMismatch as e:
                                    print(AlleleMismatch, e.args)
                                except (KeyError, ParsingError) as e:
                                    print(raw_alleles)
                                    print(variant_id)
                                    print(transcript)
                                    show_context(transcript.sequence)
                                    print(variant)

                            else:
                                print('AlleleMismatch:', raw_alleles, variant.snp_id, e.args, transcript.strand)
                                show_context(transcript.sequence)
                                continue
                        except (KeyError, ParsingError) as e:
                            print(raw_alleles)
                            print(variant_id)
                            print(transcript)
                            show_context(transcript.sequence)
                            print(variant)

                        saved = save_if_poly_a_related(by_id, transcript, alleles, variant_id)
                        #by_id[variant_id] = variant

                    if not saved and not (variant.affected_transcripts - checked_transcripts):
                        del by_id[variant_id]

                    # TODO: sanity check
                    # assert v.chr_strand == transcript.strand

    work = manager.Queue(num_workers)
    pool = []
    for i in xrange(num_workers):
        p = Process(target=do_work, args=(work, by_id, to_check_later))
        p.start()
        pool.append(p)

    chunk_size = 100000
    with fast_gzip_read(filename) as f:
        iterator = grouper(f, chunk_size)
        for i in tqdm(iterator, total=count_lines(filename) / chunk_size):
            work.put(i)

    print('Work distributed')

    for i in xrange(num_workers):
        work.put(None)

    for p in pool:
        p.join()

    print('Joined all')
