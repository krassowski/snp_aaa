import tabix

from analyses import reporter
from analyses.gtex import GTEX_DATABASE
from analyses.spidex import spidex_get_variant
from expression_database import ExpressionDatabase
from snp_parser import SPIDEX_LOCATION
from variant import Variant
from scipy.stats import pearsonr


@reporter
def gtex_on_spidex(_):
    """
    In spidex there are only SNPs (single!)

    Definitions for GTEx (from http://www.gtexportal.org/home/documentationPage):
        The effect size of the eQTLs is defined as the slope of the linear regression,
        and is computed as the effect of the alternative allele (ALT) relative to the
        reference allele (REF) in the human genome reference GRCh37/hg19
        (i.e., the eQTL effect allele is the ALT allele).


    Definitions for SPIDEX (more in spidex/README):
        dpsi_max_tissue: The delta PSI. This is the predicted change in
                         percent-inclusion due to the variant, reported
                         as the maximum across tissues (in percent).
        dpsi_zscore: This is the z-score of dpsi_max_tissue relative to the
                     distribution of dPSI that are due to common SNP.


    korelacja pearsona
    wez wszystkie punkty
    all spidex records
    """

    effect_sizes = []
    z_scores = []

    bdb = ExpressionDatabase(GTEX_DATABASE)
    tb = tabix.open(SPIDEX_LOCATION)

    for mutation_code, tissues_and_slopes in bdb.items():
        chrom, pos, ref, alt, _ = mutation_code

        # In spidex there are only SNPs (single!)
        if len(ref) != 1 or len(alt) != 1:
            continue

        pos = int(pos)

        variant = Variant(chr_name=chrom, chrom_start=pos, chrom_end=pos)

        if tissues_and_slopes:
            records = spidex_get_variant(tb, variant)
            record = None
            for record in records:
                if record.mut_allele == alt:
                    break

            if record and record.location == 'exonic':

                for tissue, slope in tissues_and_slopes:
                    effect_sizes.append(slope)
                    z_scores.append(record.dpsi_zscore)

                # print(chrom, pos, ref, alt)
                # print(record)

        pearsonr(effect_sizes, z_scores)
