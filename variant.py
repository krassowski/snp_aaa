from output_formatter import formatter


class SlottedObject(object):

    __slots__ = ()
    # volatile attributes will be ignored during objects comparison
    __volatile_attributes__ = []

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
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
        elements = []
        for attr_name in self.__slots__:
            if attr_name not in self.__volatile_attributes__:
                element = getattr(self, attr_name, None)
                if type(element) is dict:
                    element = frozenset(element.items())
                if type(element) is list:
                    element = tuple(element)
                elements.append(element)
        return tuple(elements).__hash__()

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
        'cds_start',
        'cds_end',
        'poly_aaa',
        'sequence',
        'expression'
    )

    __volatile_attributes__ = ['poly_aaa']


class Variant(SlottedObject):

    __volatile_attributes__ = (
        'affected_transcripts',
    )

    # sequence offset
    offset = 20

    __slots__ = (
        'snp_id',
        'source',
        'chr_name',
        'chr_start',
        'chr_end',
        'chr_strand',
        'ref',
        'gene',
        'correct',
        'affected_transcripts',
        #'refseq_transcript'
    )

    def __init__(self, **kwargs):

        for attr in self.__slots__:
            setattr(self, attr, None)

        self.affected_transcripts = []

        super(Variant, self).__init__(**kwargs)

    @property
    def alts(self):
        alts = set()
        for transcript in self.affected_transcripts:
            alts.update(transcript.poly_aaa.keys())
        return alts

    def as_hgvs(self):
        for alt in self.alts:
            if self.ref == alt:
                positions = self.chr_start
                event = '{ref}>{alt}'.format(ref=self.ref, alt=alt)
            elif self.ref > alt:
                positions = '{start}_{end}'.format(start=self.chr_start, end=self.chr_start + len(self.ref))
                event = 'del'
            elif self.ref < alt:
                positions = '{start}_{end}'.format(start=self.chr_start, end=self.chr_start + len(alt))
                event = 'ins{inserted}'.format(inserted=alt)
            else:
                assert False

            return '{chrom}:g.{positions}{event}'.format(chrom=self.chr_name, positions=positions, event=event)

    def is_deletion(self, alt=None):
        ref_len = len(self.ref)
        if alt:
            return ref_len > len(alt)
        else:
            return all([ref_len > len(alt) for alt in self.alts])

    def is_insertion(self, alt=None):
        ref_len = len(self.ref)
        if alt:
            return ref_len < len(alt)
        else:
            return all([ref_len < len(alt) for alt in self.alts])

    def padded_coords(self, transcript):

        assert transcript in self.affected_transcripts

        for alt in transcript.poly_aaa.keys():

            pos = int(self.chr_start)
            ref = self.ref

            if len(alt) != len(ref):
                # right padding
                if pos == 1:
                    next_aa = transcript.sequence[self.offset + 1]
                    ref = ref + next_aa
                    alt = alt + next_aa
                # left padding
                else:
                    pos -= 1
                    prev_aa = transcript.sequence[self.offset - 1]
                    ref = prev_aa + ref
                    alt = prev_aa + alt

            yield (self.chr_name, pos, ref, alt)

    @property
    def length(self):
        return self.chr_end - self.chr_start


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
