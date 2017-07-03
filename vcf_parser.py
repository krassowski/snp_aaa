from __future__ import print_function
import vcf
from parse_variants import complement


class UnknownChromosome(Exception):
    pass


class ParsingError(Exception):

    @property
    def message(self):
        return ', '.join(self.args)


def str_or_empty(x):
    if x is None:
        return '.'
    else:
        return str(x)


class VariantCallFormatParser(object):

    def __init__(self, vcf_sources):
        """
        vcf_locations: mappings name -> location for VCF files to use
        default_source: name of default VCF file to be used
        """
        self.vcf_sources = vcf_sources
        self.readers = {}
        for source, data in vcf_sources.items():
            if not data['is_alias']:
                self.readers[source] = vcf.Reader(filename=data['path'], compressed=True)
        for source, data in vcf_sources.items():
            if data['is_alias']:
                target = data['aliased_vcf']
                vcf_sources[source] = vcf_sources[target]
                self.readers[source] = self.readers[target]

    def get_by_variant(self, variant):
        """Get vcf data for given variant from one of available VCF files.

        If source is given, the VCF file matching this source will be used;
        otherwise the source will be deduced from variant.source attr

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

        pos = [variant.chr_name, variant.chr_start - 1, variant.chr_end]

        # quite common this will happen when variant is an insertion
        if pos[1] == pos[2]:
            pos[1] -= 1

        source = variant.source

        if source in ('HGMD-PUBLIC', 'PhenCode'):
            # compare:
            # 1	2160487	CD129808	CG	C	.	.	HGMD-PUBLIC_20162;TSA=deletion;E_Phenotype_or_Disease;AA=G
            # 11	64521035	CI076992	G	G.	.	.	HGMD-PUBLIC_20162;TSA=insertion;E_Phenotype_or_Disease
            # 11	64522269	HI080013	A	A.	.	.	HGMD-PUBLIC_20162;TSA=insertion;E_Phenotype_or_Disease
            # 159993316	27511	2160488	2160488	1	153979486	HGMD_MUTATION	CD129808	1	\N	8	coding_sequence_variant,upstream_gene_variant,regulatory_region_variant	23,29	12	0	\N	\N	\N	\N	418	\N	1
            # 159984197	27504	64521037	64521036	1	153987030	HGMD_MUTATION	CI076992	1	\N	8	coding_sequence_variant,downstream_gene_variant,upstream_gene_variant	23,29	10	0	\N	\N	\N	\N	418	\N	1
            # 159984221	27504	64522271	64522270	1	154081899	HGMD_MUTATION	HI080013	1	\N	8	coding_sequence_variant,upstream_gene_variant,regulatory_region_variant	23,27,29	10	0	\N	\N	\N	\N	418	\N	1

            # 17	7579583	TP53_g.11333_11334ins1	A	A.	.	.	PhenCode_20140430;TSA=sequence_alteration
            # 160029773	27509	7579585	7579584	-1	154102753	-/PhenCode_variation	TP53_g.11333_11334ins1	1	\N	6	coding_sequence_variant,non_coding_transcript_variant,non_coding_transcript_exon_variant,5_prime_UTR_variant,intron_variant,upstream_gene_variant	23,27	18	0	\N	\N	\N	\N	\N	\N	1
            pos[1] -= 2

        return self.get_by_id(source, pos, variant.snp_id)

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

    def parse(self, vcf_data, variant_source, strand, convert_to_strand=True):
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

        All alleles are strings.
            If an allele represents deletion/insertion, appropriate allele, is given as '-'.
            A lacking allele (i.e. '.' in VCF file') is indicated by '.'.

        Returns:
            tuple: (position, reference allele, alternative alleles)
        """
        if not vcf_data:
            raise ParsingError('Lack of VCF data.')

        if len(vcf_data.ALT) == 0:
            raise ParsingError('Lack of ALT.')

        ref = str(vcf_data.REF)

        alts = list(map(str_or_empty, vcf_data.ALT))

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
                    raise ParsingError(
                        'More than one strand specified: ' + ', '.join(map(str, strands)),
                        vcf_data, strands, variant_source
                    )

                vcf_strand = -1 if strands[0] == '-' else 1
                # print(vcf_strand, strand)

                if strand != vcf_strand:
                    raise ParsingError(
                        'Given strand does not match VCF strand for',
                        vcf_data, vcf_strand, strand, variant_source
                    )

            if convert_to_strand and strand == -1:
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
