class Variant(object):

    # sequence offset
    offset = 20

    attributes = (
        'refsnp_id',
        'refsnp_source',
        'chr_name',
        'chrom_start',
        'chrom_end',
        # 'allele',  # Variant Alleles
        'allele_1',  # Ancestral allele - the most frequent allele
        'minor_allele', # the second most frequent allele
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

    __slots__ = attributes + ('ref', 'gene', 'sequence', 'alts', 'poly_aaa', 'correct', '__dict__')

    def __init__(self, args):

        args = args.split('\t')

        for i, attr in enumerate(self.attributes):
            setattr(self, attr, args[i])

    def __repr__(self):

        representation = 'Variant:'

        for attr in self.attributes:
            representation += '\n\t %s: %s' % (attr, getattr(self, attr))

        return representation

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
