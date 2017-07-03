from collections import defaultdict

from jit import jit
from parse_variants import OFFSET
from poly_a import poly_a
from variant import PolyAAAData


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


def get_poly_aaa(transcript, parsed_alleles):
    """Return poly_aaa data for AffectedTranscript it will be poly_aaa related"""
    poly_aaa = analyze_poly_aaa(transcript.sequence, parsed_alleles[1:])

    for poly_a in poly_aaa.values():
        if poly_a.has or poly_a.will_have:
            print('Poly(A) detected in:', transcript.sequence)
            return poly_aaa


def show_context(seq, start=OFFSET, end=-OFFSET):
    print(seq[:start] + '>' + seq[start:end] + '<' + seq[end:])