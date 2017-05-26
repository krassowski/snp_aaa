from output_formatter import formatter


class SlottedObject(object):

    __slots__ = ()
    # volatile attributes will be ignored during objects comparison
    __volatile_attributes__ = []

    def __init__(self, **kwargs):
        for key, value in kwargs.iteritems():
            setattr(self, key, value)

    def __eq__(self, other):
        return (
            isinstance(other, self.__class__) and
            self.__slots__ == other.__slots__ and
            all(
                [
                    getattr(self, attr_name, None) == getattr(other, attr_name, None)
                    for attr_name in self.__slots__
                    if attr_name not in self.__volatile_attributes__
                ]
            )
        )

    def __hash__(self):
        return tuple(
            getattr(self, attr_name, None)
            for attr_name in self.__slots__
            if attr_name not in self.__volatile_attributes__
        ).__hash__()

    def __repr__(self):
        kwargs = {}

        for attr in self.__slots__:
            value = getattr(self, attr, None)
            if value is not None:
                kwargs[attr] = value

        return formatter(kwargs, name=self.__class__.__name__)

    def __getstate__(self):
        return [(k, getattr(self, k, None)) for k in self.__slots__]

    def __setstate__(self, data):
        for k, v in data:
            setattr(self, k, v)


class AffectedTranscript(SlottedObject):

    __slots__ = (
        'ensembl_id',
        'strand',
        'cdna_start',
        'cdna_end',
        'cds_start',
        'cds_end',
        'poly_aaa'
    )

    __volatile_attributes__ = ['poly_aaa']



class BiomartVariant(SlottedObject):

    attributes = (
        'refsnp_id',
        'refsnp_source',
        'chr_name',
        'chrom_start',
        'chrom_end',
        # 'allele',  # Variant Alleles
        'allele_1',  # Ancestral allele - the most frequent allele
        # 'minor_allele',     # the second most frequent allele
        'chrom_strand',
        'cdna_start',
        'cdna_end',
        'ensembl_gene_stable_id',
        'ensembl_transcript_stable_id',
        'ensembl_transcript_chrom_strand',
        'cds_start',
        'cds_end'
        # 'consequence_type_tv',
        # 'consequence_allele_string'
    )

    __slots__ = attributes

    def __init__(self, line_args=None, **kwargs):

        if line_args:
            args = line_args.split('\t')
            kwargs = dict(zip(self.attributes, args))

        super(BiomartVariant, self).__init__(**kwargs)

    def extract_transcript(self):
        return AffectedTranscript(
            ensembl_id=self.ensembl_transcript_stable_id,
            strand=self.ensembl_transcript_chrom_strand,
            cdna_start=self.cdna_start,
            cdna_end=self.cdna_end,
            cds_start=self.cds_start,
            cds_end=self.cds_end
        )


class Variant(SlottedObject):

    __volatile_attributes__ = (
        'affected_transcripts',
    )

    # sequence offset
    offset = 20

    biomart_attributes = (
        'refsnp_id',
        'refsnp_source',
        'chr_name',
        'chrom_start',
        'chrom_end',
        # 'allele',  # Variant Alleles
        'allele_1',  # Ancestral allele - the most frequent allele
        'minor_allele',     # the second most frequent allele
        'chrom_strand',
        'ensembl_gene_stable_id',
        # 'consequence_type_tv',
        # 'consequence_allele_string'
    )

    attributes = (
        'ref',
        'gene',
        # 'sequence',
        'sequences',
        'alts',
        'correct',
        'affected_transcripts',
        'refseq_transcript'
    )

    __slots__ = attributes + biomart_attributes

    def __init__(self, **kwargs):

        for attr in self.attributes:
            setattr(self, attr, None)

        self.affected_transcripts = set()

        super(Variant, self).__init__(**kwargs)

    def as_hgvs(self):
        # TODO: what to do for multiple alts?
        for alt in self.alts:
            if self.ref == alt:
                positions = self.chrom_start
                event = '{ref}>{alt}'.format(ref=self.ref, alt=alt)
            elif self.ref > alt:
                positions = '{start}_{end}'.format(start=self.chrom_start, end=self.chrom_start + len(self.ref))
                event = 'del'
            elif self.ref < alt:
                positions = '{start}_{end}'.format(start=self.chrom_start, end=self.chrom_start + len(alt))
                event = 'ins{inserted}'.format(inserted=alt)
            else:
                assert False

            return '{chrom}:g.{positions}{event}'.format(chrom=self.chr_name, positions=positions, event=event)

    #@property
    def is_deletion(self, alt=None):
        ref_len = len(self.ref)
        if alt:
            return ref_len > len(alt)
        else:
            return all([ref_len > len(alt) for alt in self.alts])

    #@property
    def is_insertion(self, alt=None):
        ref_len = len(self.ref)
        if alt:
            return ref_len < len(alt)
        else:
            return all([ref_len < len(alt) for alt in self.alts])

    @property
    def padded_coords(self):

        for alt in self.alts:

            pos = int(self.chrom_start)
            ref = self.ref

            if len(alt) != len(ref):
                # right padding
                if pos == 1:
                    next_aa = self.sequence[self.offset + 1]
                    ref = ref + next_aa
                    alt = alt + next_aa
                # left padding
                else:
                    pos -= 1
                    prev_aa = self.sequence[self.offset - 1]
                    ref = prev_aa + ref
                    alt = prev_aa + alt

            yield (self.chr_name, pos, ref, alt)

    @property
    def length(self):
        return self.chrom_end - self.chrom_start


class PolyAAAData(object):

    __slots__ = ('has', 'will_have', 'before', 'after')

    @property
    def change(self):
        return self.after - self.before

    @property
    def increased(self):
        return self.change > 0

    @property
    def decreased(self):
        return self.change < 0
