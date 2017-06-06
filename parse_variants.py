# -*- coding: utf-8 -*-
from __future__ import print_function
import signal
from collections import defaultdict

from poly_a import poly_a
from variant import PolyAAAData


OFFSET = 20


class ParsingError(Exception):
    pass


def decode_hgvs_code(code):
    """
    Comply to HGVS recommendations: http://www.hgvs.org/mutnomen/recs.html

    tests = [
        'FANCD1:c.658_659delGT',
        'FANCD1:c.5682C>G',
        'FANCD1:c.6275_6276delTT',
        'FANCD1:c.8504C>A',
        'FANCD1:c.5609_5610delTCinsAG',
        'FANCD1:c.1813dupA',
        'FANCD1:c.8219T>A',
        'FANCD1:c.9672dupA'
    ]

    Returns:
        gene, pos, ref, alt
    """
    # TODO testy
    ref, alt = '', ''
    gene, location = code.split(':')
    pos_type = location[0]

    match = re.match(
        '([\d]+)(_[\d]+)?([ACTG]+)?(dup|>|del|ins)([ACTG]+)(dup|>|del|ins)?([ACTG]+)?',
        location[2:]
    )

    if not match:
        raise ValueError('Cannot understand mutation code: %s' % code)

    # genomic and mitochondrial positions can be validated easily
    if pos_type in ('g', 'm'):
        pos = match.group(1)
    elif pos_type in ('c', 'n', 'p'):
        # TODO
        pos = int(match.group(1))
    else:
        raise ParsingError(
            'Wrong type of variant position specification: %s' % pos_type
        )

    event = match.group(4)
    second_event = match.group(6)

    if second_event:
        if event == 'del' and second_event == 'ins':
            ref = match.group(5)
            alt = match.group(7)
        else:
            raise ParsingError(
                'Unknown complicated event: %s and then %s' % (
                    event,
                    second_event
                )
            )
    else:
        if event == '>':
            ref = match.group(3)
            alt = match.group(5)
        elif event in ('dup', 'ins'):
            alt = match.group(5)
        elif event == 'del':
            ref = match.group(5)

    return gene, pos, ref, alt


def determine_mutation(variant, vcf_parser, offset):
    if variant.source == 'PhenCode':
        print('PhenCode variant:')
        print(variant)
        alt_source = 'PhenCode'
        gene, pos, ref, alt = decode_hgvs_code(variant.snp_id)
        alts = [alt]
    else:
        alt_source = 'VCF'

        try:
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
            strand = list(variant.affected_transcripts)[0].strand
            for transcript in variant.affected_transcripts:
                if strand != transcript.strand:
                    print('Stand mismatch for')
                    print(variant)
            pos, ref, alts = vcf_parser.parse(analysed_vcf_record, variant.source, strand)
            gene = vcf_parser.get_gene(analysed_vcf_record)
        except ParsingError as e:
            print(
                'Skipping variant: %s from %s:' % (
                    variant.snp_id,
                    variant.source
                ),
                end=' '
            )
            print(e.message)
            variant.correct = False
            return False

    # check ref sequence
    if variant.sequences:

        pos_correct = 0

        for transcript, sequence in variant.sequences.iteritems():
            seq_ref = sequence[offset:-offset]

            if ref != seq_ref:

                # Cosmic represents insertions AND bigger deletions as a range:
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=111380 - deletion (start = end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=1223920 - substitution (start = end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=4699526 - insertion (start != end)
                #   http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=4612790 - big deletion (start != end)

                if variant.source == 'COSMIC' and transcript.cds_end != transcript.cds_start:
                    # as all deletions pass ref == seq_ref test, and insertions fall there,
                    # it is needed to assure that those are corrected here. Firstly, make sure
                    # that it is an insertion:
                    assert all(len(alt) > len(ref) for alt in alts)

                    del variant.sequences[transcript]

                    # sequence was too big by 'diff'
                    diff = transcript.cds_end - transcript.cds_start + 1
                    if not pos_correct:
                        pos_correct = diff - 1
                    else:
                        assert pos_correct == diff - 1

                    sequence = sequence[:-diff]
                    transcript.cds_end = transcript.cds_start
                    assert sequence[offset:-offset] == ref

                    variant.sequences[transcript] = sequence

                else:

                    print(
                        '%s says ref is %s, but sequence analysis pointed to %s for %s, %s'
                        % (alt_source, ref, seq_ref, variant.snp_id, transcript.ensembl_id)
                    )
        pos -= pos_correct

    return {
        'gene': gene,
        'ref': ref,
        'chr_start': pos,
        'alts': alts
    }


def get_poly_a(ref_seq, alts, offset):

    has_aaa, before_len = poly_a(
        ref_seq,
        offset,
        len(ref_seq) - offset
    )

    poly_aaa = defaultdict(PolyAAAData)

    for alt in alts:

        mutated_seq = ref_seq[:offset] + str(alt) + ref_seq[-offset:]

        will_have, after_len = poly_a(
            mutated_seq,
            offset,
            len(mutated_seq) - offset
        )

        poly_aaa[alt].has = has_aaa
        poly_aaa[alt].will_have = will_have
        poly_aaa[alt].before = before_len
        poly_aaa[alt].after = after_len

    return poly_aaa


def init_worker():
    """Initialization with this function allows to terminate app

    with Ctrl+C even when multiprocessing computation is ongoing.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def get_unique_variants(variants):
    """Get only unique variants, with respect to:
        - position in genome,
        - set of alternative alleles
        - ensembl gene id

    All transcript affected by variants sharing listed properties
    will expand "affected_transcripts" set of returned variant.
    """
    # adding and removing None is workaround for https://github.com/numba/numba/issues/2152
    unique_variants = set()

    for variant in variants:
        if variant not in unique_variants:
            unique_variants.add(variant)
        else:

            for old in unique_variants:
                if old == variant:
                    old.affected_transcripts.update(variant.affected_transcripts)
                    break

    #unique_variants.remove(None)

    return unique_variants


def complement(seq):
    basic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '.': '.', '-': '-'}
    try:
        return ''.join([basic[n] for n in seq])
    except KeyError:
        # http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
        IUPAC = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A', 'Y': 'R',
                 'R': 'Y', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
        return ''.join([IUPAC[n] for n in seq])
