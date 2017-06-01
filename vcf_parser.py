from __future__ import print_function
import vcf
from berkley_hash_set import BerkleyHashSet
from fasta_sequence_db import complement


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

    def __init__(self, vcf_sources, default_source=None):
        """
        vcf_locations: mappings name -> location for VCF files to use
        default_source: name of default VCF file to be used
        """
        self.default_source = default_source
        self.vcf_sources = vcf_sources
        self.readers = {}
        for source, data in vcf_sources.iteritems():
            if not data['is_alias']:
                self.readers[source] = vcf.Reader(filename=data['path'])
        for source, data in vcf_sources.iteritems():
            if data['is_alias']:
                self.readers[source] = self.readers[data['aliased_vcf']]

    def get_by_transcript(self, transcript, pos, source):

        hgnc_name = get_hgnc(transcript.ensembl_id)

        record_id = ''.join([hgnc_name, ':c.', transcript.cds_start])

        data = []
        for record in self.get_by_pos(source, pos):
            if record.ID.startswith(record_id):
                data.append(record)

        return data

    def get_by_variant(self, variant, source=None, require_id_match=True):
        """Get vcf data for given variant from one of available VCF files.

        If source is given, the VCF file matching this source will be used;
        otherwise the source will be deduced from variant.refsnp_source attr

        Returns:
            list of vcf records as returned by vcf.fetch()
        """
        # vcf.parser uses 0-based coordinates:
        # http://pyvcf.readthedocs.org/en/latest/_modules/vcf/parser.html?highlight=coordinates
        # ensembl uses 1-based coordinates:
        # http://www.ensembl.org/info/docs/api/core/core_tutorial.html#coordinates

        # (for performance reasons incorporated into a list creation, left as a comment)
        # to get to vcf stored data by vcf reader, change coordinates to 0-based
        # pos[1] -= 1
        # pos[2] -= 1
        # and represent them as a range
        # pos[2] += 1

        pos = [variant.chr_name, variant.chrom_start - 1, variant.chrom_end]

        # quite common this will happen when variant is an insertion
        if pos[1] == pos[2]:
            pos[1] -= 1

        if not source:
            source = variant.refsnp_source

        if source not in self.readers:
            print('Unknown source: %s' % source)
            print(variant)

            source = self.default_source

            if not source:
                raise ParsingError(
                    'Either source or default_source must be present'
                    ' to choose vcf file to look up'
                )

        if require_id_match:
            return self.get_by_id(source, pos, variant.refsnp_id)
        else:
            return self.get_by_pos(source, pos)

    def get_by_pos(self, reader_name, pos):
        """Ultimately all record getters end up reaching for get_by_pos
        as this is the way the reader retrieve data from file using tabix.
        """
        try:
            reader = self.readers[reader_name]
        except KeyError:
            raise ParsingError(
                '%s is not a known source for VCF file. '
                'Known sources are: %s.' % (
                    reader_name,
                    self.readers.keys()
                )
            )

        vcf_data = []

        try:
            for record in reader.fetch(*pos):
                vcf_data.append(record)
        except ValueError:
            raise ParsingError('Failed to create iterator for pos: %s' % pos)

        return vcf_data

    def get_by_id(self, reader_name, pos, record_id):

        vcf_data = []

        for record in self.get_by_pos(reader_name, pos):
            if record.ID == record_id:
                vcf_data.append(record)

        return vcf_data

    def parse(self, vcf_data, variant_source, strand):
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
            tuple: (position, reference allele, alternative alleles)
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

        if self.vcf_sources[variant_source]['given_as_positive_strand_only']:
            # take a look at this example:
            # 1	91404597	COSM913148	AAGAATTT	A	.	.	GENE=ZNF644;STRAND=-;CDS=c.2307_2313delAAATTCT;AA=p.L769fs*36;CNT=1
            # first left padding is removed (before this if) (so we have AGAATTT),
            # then we get reversed complement to obtain AAATTCT
            if 'STRAND' in vcf_data.INFO:
                strands = vcf_data.INFO['STRAND']

                if len(strands) > 1:
                    print('More than one strand specified: ' + ', '.join(map(str, strands)))
                    print(vcf_data, strands, variant_source)

                vcf_strand = -1 if strands[0] == '-' else 1

                if strand != vcf_strand:
                    print('Given strand does not match VCF strand for')
                    print(vcf_data, vcf_strand, strand, variant_source)

            if strand == -1:
                ref = complement(ref)[::-1]
                alts = [complement(alt)[::-1] for alt in alts]

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
