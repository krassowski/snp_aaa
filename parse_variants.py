# -*- coding: utf-8 -*-
from __future__ import print_function
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
