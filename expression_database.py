import os

from tqdm import tqdm

from analyses.spidex import convert_to_strand
from berkley_hash_mutable_callback import BerkleyHashSet, BerkleyHashList
from multiprocess import fast_gzip_read, count_lines

DEFAULT_PATH = 'GTEx_Analysis_v6p_eQTL'
DEFAULT_SUFFIX = '_Analysis.v6p.signif_snpgene_pairs.txt.gz'
DEFAULT_GENE_SUFFIX = '_Analysis.v6p.egenes.txt.gz'


GTEX_TISSUES = [
    'Adipose_Subcutaneous',
    'Adipose_Visceral_Omentum',
    'Adrenal_Gland',
    'Artery_Aorta',
    'Artery_Coronary',
    'Artery_Tibial',
    'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia',
    'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum',
    'Brain_Cortex',
    'Brain_Frontal_Cortex_BA9',
    'Brain_Hippocampus',
    'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia',
    'Brain_Putamen_basal_ganglia',
    'Breast_Mammary_Tissue',
    'Cells_EBV-transformed_lymphocytes',
    'Cells_Transformed_fibroblasts',
    'Colon_Sigmoid',
    'Colon_Transverse',
    'Esophagus_Gastroesophageal_Junction',
    'Esophagus_Mucosa',
    'Esophagus_Muscularis',
    'Heart_Atrial_Appendage',
    'Heart_Left_Ventricle',
    'Liver',
    'Lung',
    'Muscle_Skeletal',
    'Nerve_Tibial',
    'Ovary',
    'Pancreas',
    'Pituitary',
    'Prostate',
    'Skin_Not_Sun_Exposed_Suprapubic',
    'Skin_Sun_Exposed_Lower_leg',
    'Small_Intestine_Terminal_Ileum',
    'Spleen',
    'Stomach',
    'Testis',
    'Thyroid',
    'Uterus',
    'Vagina',
    'Whole_Blood'
]


class ExpressedGenes(BerkleyHashList):
    pass


class ExpressionDatabase(BerkleyHashSet):

    def get_by_mutation(self, mutation, transcript):

        results = {}

        for chrom, pos, ref, alt in mutation.padded_coords(transcript):

            strand = '+' if transcript.strand == 1 else '-'

            key = '_'.join(map(str, [
                chrom,
                pos,
                convert_to_strand(ref, strand),
                convert_to_strand(alt, strand)
            ])) + '_b37'

            # https://www.gtexportal.org/home/snp/rs7841915
            # if mutation.snp_id == 'rs7841915':
            #   print(key)

            data = self[key]

            results[alt] = [
                datum.split(',')
                for datum in data
            ]

        return results

    def items(self):
        for key, data in super(ExpressionDatabase, self).items():
            yield key.split('_'), [
                (t, float(e), g)
                for datum in data
                for (t, e, g) in [datum.split(',')]
            ]


def expression_file_paths(tissues_list, path, suffix):
    for tissue_name in tissues_list:
        file_name = tissue_name + suffix
        file_path = os.path.join(path, file_name)
        yield file_path


def count_all(tissues_list, path=DEFAULT_PATH, suffix=DEFAULT_SUFFIX):
    return sum(
        count_lines(path)
        for path in expression_file_paths(tissues_list, path, suffix)
    )


def import_expressed_genes(
        bdb,
        tissues=GTEX_TISSUES,
        path=DEFAULT_PATH,
        suffix=DEFAULT_GENE_SUFFIX
):
    print('Importing expressed genes:')

    count = count_all(tissues, path, suffix)

    with tqdm(total=count) as progress:

        for tissue_name in tissues:
            file_name = tissue_name + suffix
            file_path = os.path.join(path, file_name)
            print('Loading', file_name)

            with fast_gzip_read(file_path) as file_object:
                # skip header
                next(file_object)

                for line in file_object:
                    data = line.split()
                    """
                    gene_id: (
                        gene_name,
                        gene_chr,
                        gene_start,
                        gene_end,
                        strand,
                    )
                    """
                    if not bdb[data[0]]:
                        bdb[data[0]].extend(data[1:6])
                    else:
                        assert bdb[data[0]] == data[1:6]

                    progress.update(1)


#@jit
def iterate_over_expression(
    tissues_list=GTEX_TISSUES,
    path=DEFAULT_PATH,
    suffix=DEFAULT_SUFFIX
):
    for tissue_name in tissues_list:
        file_name = tissue_name + suffix
        file_path = os.path.join(path, file_name)
        print('Loading', file_name)

        with fast_gzip_read(file_path, processes='all') as file_object:

            header_line = next(file_object)
            header = dict()

            for position, name in enumerate(header_line.split()):
                header[name] = position

            slope_pos = header['slope']
            gene_id_pos = header['gene_id']
            variant_id_pos = header['variant_id']

            for line in file_object:
                data = line.split()
                yield (data[variant_id_pos], tissue_name, data[slope_pos], data[gene_id_pos])


def import_expression_data(
    bdb,
    tissues=GTEX_TISSUES,
    path=DEFAULT_PATH,
    suffix=DEFAULT_SUFFIX
):
    """Import expression data from GTEx portal files "signif_snpgene_pairs"
    into given database

    Args:
        bdb: database object with >dict of sets< interface (ExpressionDatabase)
        path: path to folder with signif_snpgene_pairs.txt.gz files
    """
    print('Importing expression data:')

    count = count_all(tissues, path, suffix)

    expression_data = iterate_over_expression(tissues, path, suffix)

    for expr in tqdm(expression_data, total=count):
        bdb[expr[0]].add(','.join(expr[1:]))

