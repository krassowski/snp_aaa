from pyfaidx import Fasta

import multiprocess
from jit import jit
from multiprocess import get_manager
from parse_variants import OFFSET
from settings import TRANSCRIPT_DB_PATH
from variant import AffectedTranscript
from variant_sources.ensembl.poly_aaa import get_poly_aaa, show_context
from variant_sources.ensembl.allele_parsing import AlleleMismatch, parse_alleles
from vcf_parser import ParsingError

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


def do_work(progress, in_queue, static_args, accepted_pairs, to_check_pairs):

    path, transcript_strand, accepted_consequences, limit_to_transcripts, only_poly_a = static_args

    transcripts_db = Fasta(TRANSCRIPT_DB_PATH, key_function=take_transcript_id_without_version)

    @jit
    def good_consequence(data):
        # to restrict to somatic mutations, use: and data[somatic_pos] == '1'
        return not accepted_consequences.isdisjoint(data[consequence_pos].split(','))

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


def load_poly_a_transcript_variant_pairs(path, transcript_strand, accepted_consequences, limit_to_transcripts, only_poly_a):

    manager = get_manager()
    accepted_transcript_variant_pairs = manager.list()
    to_check_transcript_variant_pairs = manager.list()

    filename = path + 'transcript_variation.txt.gz'

    multiprocess.parse_gz_file(
        filename,
        do_work,
        static_args=[(path, transcript_strand, accepted_consequences, limit_to_transcripts, only_poly_a)],
        shared_args=[accepted_transcript_variant_pairs, to_check_transcript_variant_pairs],
        chunk_size=1000000
    )

    return accepted_transcript_variant_pairs, to_check_transcript_variant_pairs


@jit
def take_transcript_id_without_version(full_id):
    """Returns transcript id without version and everything which is after version separating comma.

    Example:
        Input: ESNT_0001.4
        Output: ESNT_0001

        Input: ESNT_0002.2.some_annotation
        Output: ESNT_0002

    """
    return full_id.split('.')[0]


@jit
def int_or_none(x):
    if x == '\\N':
        return None
    return int(x)
