import tabix

from analyses import reporter
from analyses.spidex import spidex_get_variant, choose_record
from settings import SPIDEX_LOCATION
from variant import Variant


@reporter
def zcrb1(_):
    """
    Data based on:
    http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=109189

    In spidex ref and alt is always from forward strand: (!)
        ref_allele: The reference allele at the variant position (forward-strand)
        mut_allele: The mutant allele at the variant position (forward-strand)
    """
    ref = 'G'
    alt = 'A'
    chrom = '12'
    strand = '-'
    pos = 42707711

    tb = tabix.open(SPIDEX_LOCATION)

    variant = Variant(
        chr_name=chrom,
        chr_start=pos,
        chr_end=pos,
        chr_strand=strand,
        snp_id='ZCRB1:c.411G>A',
        ref=ref
    )

    records = spidex_get_variant(tb, variant)
    record = choose_record(records, variant, alt, location='exonic', strict=True)

    assert record

    z_score = float(record.dpsi_zscore)

    print(z_score, record)

