from __future__ import print_function

import gc
from collections import defaultdict

from pyfaidx import Fasta
from tqdm import tqdm

import multiprocess
from commands import SourceSubparser
from multiprocess import fast_gzip_read, count_lines, manager
from parse_variants import OFFSET, complement
from poly_a import poly_a
from snp_parser import jit
from snp_parser import vcf_mutation_sources, TRANSCRIPT_DB_PATH
from variant import Variant, AffectedTranscript, PolyAAAData
from variant_sources import variants_getter
from vcf_parser import VariantCallFormatParser, ParsingError


def show_context(seq, start=OFFSET, end=-OFFSET):
    print(seq[:start] + '>' + seq[start:end] + '<' + seq[end:])


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
    #default='/media/ramdisk/head/'
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
    acction='store_true',
    help=(
        'include all variants rather than only those which are poly_a related'
    )
)


@jit
def analyze_poly_aaa(sequence, alternative_alleles):

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


def parse_alleles(alleles, transcript):
    sequence_ref_allele = transcript.sequence[OFFSET:-OFFSET]

    if not alleles:
        print('No alleles for:')
        print(transcript)
        return

    if len(alleles) == 2:
        allele = alleles[1]

        if allele.startswith('(') and allele.endswith(')'):
            alleles = parse_short_tandem_repeat(allele, alleles)

    alleles = ['' if a == '-' else a for a in alleles]

    if transcript.strand == -1:
        alleles = [
            complement(alt)[::-1] for alt in alleles
        ]

    if sequence_ref_allele != alleles[0]:
        raise AlleleMismatch(sequence_ref_allele, alleles[0])

    return alleles


def parse_short_tandem_repeat(allele, alleles):
    r"""
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


def get_poly_aaa(transcript, parsed_alleles):
    """Return poly_aaa data for AffectedTranscript it will be poly_aaa related"""
    poly_aaa = analyze_poly_aaa(transcript.sequence, parsed_alleles[1:])

    for poly_a in poly_aaa.values():
        if poly_a.has or poly_a.will_have:
            print('Poly(A) detected in:', transcript.sequence)
            return poly_aaa


def aggregate_transcripts(transcript_variant_pairs):
    print('Aggregating transcripts')
    transcripts_by_id = defaultdict(list)
    for internal_variant_id, transcript, alleles in tqdm(transcript_variant_pairs):
        transcripts_by_id[internal_variant_id].append((transcript, alleles))
    return transcripts_by_id


@variants_getter
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
    accepted_transcripts_by_variant_id = aggregate_transcripts(accepted_transcript_variant_pairs)
    to_check_transcripts_by_variant_id = aggregate_transcripts(to_check_transcript_variant_pairs)

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


def load_poly_a_transcript_variant_pairs(path, transcript_strand, accepted_consequences, limit_to_transcripts, only_poly_a=True):

    filename = path + 'transcript_variation.txt.gz'

    headers = [
        'transcript_variation_id', 'variation_feature_id', 'feature_stable_id',
        'allele_string',
        'somatic',
        'consequence_types',
        'cds_start', 'cds_end',
        'cdna_start', 'cdna_end',
        'translation_start', 'translation_end',
        'distance_to_transcript',
        'codon_allele_string',
        'pep_allele_string',
        'hgvs_genomic', 'hgvs_transcript', 'hgvs_protein',
        'polyphen_prediction', 'polyphen_score',
        'sift_prediction', 'sift_score',
        'display'
    ]
    # somatic_pos = headers.index('somatic')
    consequence_pos = headers.index('consequence_types')

    @jit
    def good_consequence(data):
        # to restrict to somatic mutations, use: and data[somatic_pos] == '1'
        return not accepted_consequences.isdisjoint(data[consequence_pos].split(','))

    def do_work(progress, in_queue, accepted_pairs, to_check_pairs):

        transcripts_db = Fasta(TRANSCRIPT_DB_PATH, key_function=lambda x: x.split('.')[0])

        @jit
        def get_any_sequence(transcript, raw_start, raw_end):

            if not (raw_start and raw_end):
                # print('No sequence coordinates for %s transcript' % transcript.ensembl_id)
                return

            if transcript not in transcripts_db:
                return

            if raw_end - raw_start > 2 * OFFSET:
                # print('Skipping transcript %s of variant %s: too wide variation' % (transcript.ensembl_id, variant.snp_id))
                return

            whole_seq = str(transcripts_db[transcript])

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

        if limit_to_transcripts:
            # when limiting transcripts number significantly, additional overhead of one function is not a problem
            @jit
            def get_sequence(transcript, raw_start, raw_end):
                if transcript in limit_to_transcripts:
                    return get_any_sequence(transcript, raw_start, raw_end)
        else:
            # otherwise it's better to limit overhead as much as possible
            get_sequence = get_any_sequence

        for lines in multiprocess.parser(progress, in_queue):
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
                cds_start = int_or_none(data[6])
                cds_end = int_or_none(data[7])

                sequence = get_sequence(ensembl_id, cds_start, cds_end)

                if not sequence:
                    continue

                transcript = AffectedTranscript(
                    cds_start=cds_start,
                    cds_end=cds_end,
                    strand=transcript_strand[ensembl_id],
                    ensembl_id=ensembl_id,
                    sequence=sequence
                )

                internal_variant_id = int(data[1])

                if len(alleles) < 1:
                    # 152791410	27516	66766360	66766362	1	153751119	(GGC)18	rs869109080	1	\N	1	non_coding_transcript_variant,non_coding_transcript_exon_variant,NMD_transcript_variant,inframe_insertion,regulatory_region_variant	1	16	0	\N	\N	\N	1	\N	\N	0
                    continue

                # TODO: why some PhenCode variations ends with 'puh'?
                if alleles[1] in ('COSMIC_MUTATION', 'HGMD_MUTATION', 'PhenCode_variation') or alleles[1].endswith('puh'):
                    to_check_pairs.append((internal_variant_id, transcript, alleles))
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
                    to_check_pairs.append((internal_variant_id, transcript, alleles))
                    continue
                except (KeyError, ParsingError) as e:
                    print(e)
                    print(alleles)
                    print(internal_variant_id)
                    print(transcript)
                    show_context(transcript.sequence)
                    continue

                poly_aaa = get_poly_aaa(transcript, alleles)

                if poly_aaa or not only_poly_a:
                    transcript.poly_aaa = poly_aaa
                    accepted_pairs.append((internal_variant_id, transcript, alleles))

    accepted_transcript_variant_pairs = manager.list()
    to_check_transcript_variant_pairs = manager.list()

    multiprocess.parse_gz_file(
        filename,
        do_work,
        shared_args=[accepted_transcript_variant_pairs, to_check_transcript_variant_pairs],
        chunk_size=1000000
    )

    return accepted_transcript_variant_pairs, to_check_transcript_variant_pairs


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


def load_variants_details(path, seq_region, sources, all_accepted_by_id, all_to_check_by_id, only_poly_a):
    """
    feature_table = [
        'variation_feature_id',
        'seq_region_id', 'seq_region_start', 'seq_region_end', 'seq_region_strand',
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
        'minor_allele', 'minor_allele_freq', 'minor_allele_count',
        'alignment_quality',
        'evidence_attribs',
        'clinical_significance',
        'display',
    ]
    """
    filename = path + 'variation_feature.txt.gz'

    def test_poly_a_and_amend_transcript_if_present(transcript, variant, alleles, only_poly_a=True):
        try:
            pos, ref, alts = get_alt_alleles(variant, transcript.strand, convert_to_strand=False)
        except ParsingError as e:
            if variant.source == 'COSMIC':
                # TODO: more specific exception
                # some COSMIC mutations are mapped twice to ensembl but only one
                # position is given in Cosmic VCF file. In such cases, I stick with
                # Cosmic VCF file. There is no need to print out such cases
                pass
            else:
                print('ParsingError:', e.message, e.args)
                print(variant.snp_id, alleles)
            return

        if any('.' in alt for alt in alts):
            if variant.source == 'HGMD-PUBLIC':
                pass
                # TODO: count this
                # print('Skipping variant from HGMD without complete allele data: %s' % v.snp_id)
            else:
                print('ParsingError: IncorrectAllele, with dot, not from HGMD')
            return

        raw_alleles = [ref] + alts

        try:
            alleles = parse_alleles(raw_alleles, transcript)
        except AlleleMismatch as e:
            if variant.source == 'COSMIC':

                # Cosmic represents insertions AND bigger deletions as a range (cds start != cds end):
                # though other sources represent insertions as a point event (cds start == cds end).
                # This causes COSMIC insertions to be wrongly recognised and AlleleMismatch to be raised.
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=111380 - deletion (start = end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=1223920 - substitution (start = end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=4699526 - insertion (start != end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=4612790 - big deletion (start != end)

                # Let's handle insertions:
                if ref == '' and all(len(alt) > len(ref) for alt in alts):
                    # sequence was too big by 'diff'
                    diff = transcript.cds_end - transcript.cds_start + 1

                    sequence = transcript.sequence[:-diff]
                    transcript.cds_end = transcript.cds_start
                    assert sequence[OFFSET:-OFFSET] == ref

                    transcript.sequence = sequence

                else:
                    # sometimes ensembl has wrongly stated ref allele too.
                    print('AlleleMismatch:', e.args)
                    print(variant)
                    print(transcript)
                    show_context(transcript.sequence)
                    return

                try:
                    alleles = parse_alleles(raw_alleles, transcript)
                except AlleleMismatch as e:
                    print('AlleleMismatch:', e.args)
                    return
                except (KeyError, ParsingError) as e:
                    print(raw_alleles)
                    print(transcript)
                    show_context(transcript.sequence)
                    print(variant)
                    return

            else:
                print('AlleleMismatch:', raw_alleles, variant.snp_id, e.args, transcript.strand)
                show_context(transcript.sequence)
                return

        except (KeyError, ParsingError) as e:
            print(raw_alleles)
            print(transcript)
            show_context(transcript.sequence)
            print(variant)
            return

        poly_aaa = get_poly_aaa(transcript, alleles)

        if poly_aaa or not only_poly_a:
            transcript.poly_aaa = poly_aaa
            return True

    def do_work(progress, in_queue, accepted_by_id, to_check_by_id, variants_by_id_out):
        # now I have a guarantee that a single ID will come only in this process (!)

        for lines in multiprocess.parser(progress, in_queue):

            variant_by_id_buffer = {}

            for line in lines:

                if line is None:
                    in_queue.task_done()
                    return

                data = line.split('\t')
                variant_id = int(data[0])

                in_accepted = variant_id in accepted_by_id
                in_to_check = variant_id in to_check_by_id

                if not (in_accepted or in_to_check):
                    continue

                variant = Variant(
                    chr_name=seq_region[int(data[1])],
                    chr_start=int(data[2]),
                    chr_end=int(data[3]),
                    chr_strand=data[4],
                    snp_id=data[7],
                    source=sources[int(data[10])]
                )

                if in_accepted:
                    variant.affected_transcripts = [
                        transcript_and_allele[0] for transcript_and_allele in accepted_by_id[variant_id]
                    ]

                if in_to_check:
                    newly_accepted_transcripts = []
                    for transcript, alleles in to_check_by_id[variant_id]:
                        accepted = test_poly_a_and_amend_transcript_if_present(transcript, variant, alleles, only_poly_a)
                        if accepted:
                            newly_accepted_transcripts.append(transcript)

                    variant.affected_transcripts.extend(newly_accepted_transcripts)

                if variant.affected_transcripts:
                    ref = variant.affected_transcripts[0].sequence[OFFSET:-OFFSET]
                    if any(ref != t.sequence[OFFSET:-OFFSET] for t in variant.affected_transcripts):
                        print('Mismatch between transcripts reference allele!')
                        print(variant)
                        continue
                    variant.ref = ref
                    variant_by_id_buffer[variant_id] = variant

            variants_by_id_out.update(variant_by_id_buffer)

    variants_by_id_out = manager.dict()

    multiprocess.parse_gz_file(
        filename,
        do_work,
        static_args=[all_accepted_by_id, all_to_check_by_id],
        shared_args=[variants_by_id_out],
        chunk_size=100000
    )

    return variants_by_id_out
