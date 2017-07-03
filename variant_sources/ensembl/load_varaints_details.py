import multiprocess
from multiprocess import get_manager
from parse_variants import OFFSET
from settings import vcf_mutation_sources
from variant import Variant
from variant_sources.ensembl.allele_parsing import AlleleMismatch, parse_alleles
from variant_sources.ensembl.poly_aaa import get_poly_aaa, show_context
from vcf_parser import ParsingError, VariantCallFormatParser


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


def do_work(progress, in_queue, accepted_by_id, to_check_by_id, static_args, variants_by_id_out):
    path, seq_region, sources, all_accepted_by_id, all_to_check_by_id, only_poly_a = static_args
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

    manager = get_manager()
    variants_by_id_out = manager.dict()

    multiprocess.parse_gz_file(
        filename,
        do_work,
        static_args=[all_accepted_by_id, all_to_check_by_id, (path, seq_region, sources, all_accepted_by_id, all_to_check_by_id, only_poly_a)],
        shared_args=[variants_by_id_out],
        chunk_size=100000
    )

    return variants_by_id_out


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

