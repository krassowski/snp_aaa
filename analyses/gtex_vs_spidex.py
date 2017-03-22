import tabix

from tqdm import tqdm

from analyses import reporter
from analyses.spidex import spidex_get_variant, choose_record, convert_to_strand
from expression_database import ExpressedGenes, iterate_over_expression, count_all, TISSUES_LIST, \
    import_expressed_genes
from snp_parser import SPIDEX_LOCATION
from variant import Variant
from scipy.stats import pearsonr
import numpy as np


def create_path_for_genes_db(tissues):
    """Create almost-certainly unique path for a database which will contain
    information about genes having variant-egene pairs in given tissues."""

    from hashlib import sha224

    tissues_hash_code = sha224(','.join(tissues)).hexdigest()

    return 'genes_{hash_code}.db'.format(
        hash_code=tissues_hash_code
    )


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
    print('Gtex_on_spidex')

    effect_sizes = []
    z_scores = []

    tissues_list = TISSUES_LIST
    # Use "Brain cortex" for basic tests - it's very small
    # tissues_list = ['Brain_Cortex']
    # Use "Adipose Subcutaneous" for larger tests
    # tissues_list = ['Adipose_Subcutaneous']

    path = create_path_for_genes_db(tissues_list)
    genes = ExpressedGenes(path)
    genes.reset()

    import_expressed_genes(
        genes,
        tissues_list=tissues_list
    )

    tb = tabix.open(SPIDEX_LOCATION)

    count = count_all(tissues_list)

    with open('gtex_vs_spidex.csv', 'w') as f:
        f.write(
            'chrom,pos,ref,alt,gene_name,gene_start,gene_end,gene_strand,'
            'gtex_slope,gtex_tissue,spidex_dpsi_zscore,spidex_dpsi_max_tissue\n'
        )
        for mutation_code, tissue, slope, gene in tqdm(iterate_over_expression(tissues_list), total=count):
            chrom, pos, ref, alt, _ = mutation_code.split('_')

            # In spidex there are only SNPs (single!)
            if len(ref) != 1 or len(alt) != 1:
                continue

            pos = int(pos)

            if not genes[gene]:
                print('gene %s not present in data' % gene)
                continue
            name, chrom, start, end, strand = genes[gene]

            # only genes overlapping with given mutation
            if not (int(start) <= pos <= int(end)):
                continue

            variant = Variant(
                chr_name=chrom,
                chrom_start=pos,
                chrom_end=pos,
                chrom_strand=strand,
                refsnp_id='-',
                ref=convert_to_strand(ref, strand)
            )

            records = spidex_get_variant(tb, variant)
            record = choose_record(records, variant, convert_to_strand(alt, strand), strict=False)

            if record:
                if name != record.gene:
                    print('Gene name mismatch %s %s!' % (name, record.gene))
                    continue

                f.write(','.join(map(str, [
                    chrom, pos, convert_to_strand(ref, strand),
                    convert_to_strand(alt,strand), name, start, end, strand,
                    slope, tissue, record.dpsi_zscore, record.dpsi_max_tissue
                ])) + '\n')
                effect_sizes.append(float(slope))
                #z_scores.append(float(record.dpsi_max_tissue))
                z_scores.append(float(record.dpsi_zscore))

    print('Found %s pairs gtex-spidex' % len(effect_sizes))

    with open('effect_sizes.txt', 'w') as f:
        f.writelines([str(x) + '\n' for x in effect_sizes])

    with open('z_scores.txt', 'w') as f:
        f.writelines([str(x) + '\n' for x in z_scores])

    effect_sizes = np.array(effect_sizes)
    z_scores = np.array(z_scores)

    # "The null hypothesis is that the two variables are uncorrelated."
    pearson_coef, p_value = pearsonr(effect_sizes, z_scores)

    print(pearson_coef, p_value)

    #import code
    #code.interact(local=locals())

