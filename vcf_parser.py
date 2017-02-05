from __future__ import print_function
import vcf
from berkley_hash_set import BerkleyHashSet


class UnknownChromosome(Exception):
    pass


def get_hgnc(ensembl_transcript_id):

    hgnc_by_ensembl = BerkleyHashSet('hgnc_by_ensembl.db')

    names = hgnc_by_ensembl[ensembl_transcript_id]

    # convert set to a list for indexing support
    names = [name for name in names]

    if len(names) == 0:
        raise ValueError('No HGNC for transcript: ' + ensembl_transcript_id)
    elif len(names) > 1:
        print('Multiple HGNC identifiers for transcript: ' + ensembl_transcript_id)

    return names[0]


class ParsingError(Exception):
    pass


class VariantCallFormatParser(object):

    def __init__(self, vcf_locations):

        self.readers = {
            source: vcf.Reader(filename=path)
            for source, path in vcf_locations.items()
        }

    def get_by_variant(self, variant):

        pos = [str(variant.chr_name), variant.chrom_start, variant.chrom_end]

        # to get to vcf stored data by vcf reader, change coordinates to 0-based
        pos[1] -= 1
        pos[2] -= 1

        # and represent them as a range
        pos[2] += 1

        source = variant.refsnp_source
        record_id = variant.refsnp_id

        if source not in self.readers:
            source = 'other'

            hgnc_name = get_hgnc(variant.ensembl_transcript_stable_id)

            record_id = ''.join([
                hgnc_name, ':c.', variant.cds_start,
                variant.allele_1, '>', variant.minor_allele
            ])

        return self.get_by_id(source, pos, record_id)

    def get_by_id(self, reader_name, pos, record_id):

        reader = self.readers[reader_name]
        vcf_data = None

        for record in reader.fetch(*pos):
            if record.ID == record_id:
                vcf_data = record
                break

        return vcf_data

    @staticmethod
    def parse(vcf_data):
        """Parse VCF data, retrieved with retrieve() func, with
        compliance to VCF 4.3 specification, particularly:

        "For simple insertions and deletions in which either the REF or one of
        the ALT alleles would otherwise be null/empty, the REF and ALT.
        Strings must include the base before the event (which must be reflected
        in the POS field), unless the event occurs at position 1 on the contig
        in which case it must include the base after the event; this padding base
        is not required (although it is permitted) for e.g. complex substitutions
        or other events where all alleles have at least
        one base represented in their Strings"

        Returns:
            tuple: (position, reference alelle, alternative alleles)
        """
        if not vcf_data:
            raise ParsingError('Lack of VCF data.')

        if len(vcf_data.ALT) == 0:
            raise ParsingError('Lack of ALT.')

        ref = str(vcf_data.REF)
        alts = list(map(str, vcf_data.ALT))

        # let's recognize indel mutations and remove vcf padding from alt/ref variables
        left = 0
        right = 0

        # padding should be possible to determine from any of alternative
        # allele strings, so why not to get the first one?
        alt = alts[0]
        # but, to be certain, lets check if all alts are of the same length:
        if any(len(a) != len(alt) for a in alts):
            raise ParsingError('Alts are of different lengths: ' + str(alts))

        # left side padding, most common
        if alt[0] == ref[0]:
            left += 1
        # right side padding, pos should be one
        elif alt[-1] == ref[-1]:
            if vcf_data.POS == 1:
                right += 1
            else:
                print(
                    'Unusual padding detected (left while pos: %s) '
                    'where raw ref: "%s" and analysed alt: "%s"' %
                    (vcf_data.POS, ref, alt)
                )

        pos = vcf_data.POS + left
        ref = ref[left:-right or None]
        alts = [a[left:-right or None] for a in alts]

        # validation
        if left and right:
            raise ParsingError(
                'Wrong padding detected: both left and right present' +
                'where ref: "%s" and analysed alt: "%s"' % (ref, alt)
            )

        if not (left or right) and len(ref) != len(alt):
            raise ParsingError(
                'No padding detected despite alt/ref of different lengths ' +
                'where ref: "%s" and analysed alt: "%s"' % (ref, alt)
            )

        return pos, ref, alts

    @staticmethod
    def get_gene(vcf_data):
        try:
            gene = vcf_data.INFO['GENEINFO'].split(':')[0]
        except KeyError:
            try:
                gene = vcf_data.INFO['GENE'][0]
            except KeyError:
                gene = None
        return gene
