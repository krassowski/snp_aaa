from parse_variants import OFFSET, complement
from vcf_parser import ParsingError


class AlleleMismatch(Exception):
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
            repeat_count = (max_count + min_count) // 2
        else:
            repeat_count = int(repeat_count)
        alt = ref * repeat_count
        if alleles[0] != ref:
            raise ParsingError('STR has alleles different than deduced ref', alleles, allele, ref)
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