import tabix

from tqdm import tqdm

from analyses import reporter
from analyses.gtex import GTEX_DATABASE, GTEX_GENES
from analyses.spidex import spidex_get_variant, choose_record, convert_to_strand, StrandMismatch
from expression_database import ExpressionDatabase, ExpressedGenes
from snp_parser import SPIDEX_LOCATION
from variant import Variant
from scipy.stats import pearsonr
import numpy as np


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

        ref_allele: The reference allele at the variant position (forward-strand)
        mut_allele: The mutant allele at the variant position (forward-strand)

    """

    effect_sizes = []
    z_scores = []

    bdb = ExpressionDatabase(GTEX_DATABASE)
    genes = ExpressedGenes(GTEX_GENES)
    tb = tabix.open(SPIDEX_LOCATION)

    for mutation_code, data in tqdm(bdb.items(), total=len(bdb)):
        if not data:
            continue

        chrom, pos, ref, alt, _ = mutation_code

        # In spidex there are only SNPs (single!)
        if len(ref) != 1 or len(alt) != 1:
            continue

        pos = int(pos)

        for tissue, slope, gene in data:
            name, chrom, start, end, strand = genes[gene]

            # only those cis affected genes within given mutation lies
            if not (pos >= int(start) and pos <= int(end)):
                continue

            variant = Variant(
                chr_name=chrom,
                chrom_start=pos,
                chrom_end=pos,
                chrom_strand=strand,
                refsnp_id='unknown',
                ref=convert_to_strand(ref, strand)
            )

            records = spidex_get_variant(tb, variant)
            try:
                record = choose_record(records, variant, convert_to_strand(alt, strand), silent_intronic=True)
            except StrandMismatch:
                #if name != record.gene:
                # print('!')
                print('Strand mismatch %s %s %s ' % (chrom, pos, gene))
                continue

            # TODO: GTEx yields lots of genes in single locus. To be investigated;

            if record:

                effect_sizes.append(slope)
                z_scores.append(float(record.dpsi_zscore))

            # print(chrom, pos, ref, alt)
            # print(record)

    effect_sizes = np.array(effect_sizes)
    z_scores = np.array(z_scores)

    # "The null hypothesis is that the two variables are uncorrelated."
    pearson_coef, p_value = pearsonr(effect_sizes, z_scores)

    print(pearson_coef, p_value)
