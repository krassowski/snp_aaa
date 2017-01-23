from berkley_hash_set import BerkleyHashSet
import gzip
from tqdm import tqdm
import os
from operator import itemgetter


TISSUES_LIST = [
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


class ExpressionDatabase(BerkleyHashSet):

    def get_by_mutation(self, mutation):

        key = '_'.join(map(str, [
            mutation.chr_name,
            int(mutation.chrom_start),   # use 0 based  # TODO check off by 1?
            mutation.ref,
            mutation.alt
        ])) + '_b37'

        data = self[key]

        return [
            datum.split(',')
            for datum in data
        ]

    def __getitem__(self, mutation_key):
        value = super(ExpressionDatabase, self).__getitem__(
            mutation_key
        )
        return value

    def __setitem__(self, mutation_key, value):

        return super(ExpressionDatabase, self).__setitem__(
            mutation_key, value
        )


def import_expression_data(
    bdb,
    tissues_list=TISSUES_LIST,
    path='GTEx_Analysis_v6p_eQTL',
    suffix='_Analysis.v6p.signif_snpgene_pairs.txt.gz'
):
    """Import expression data from GTEx portal files "signif_snpgene_pairs"
    into given database

    Args:
        bdb: database object with >dict of sets< interface (ExpressionDatabase)
        path: path to folder with signif_snpgene_pairs.txt.gz files
    """
    print('Importing expression data:')

    for tissue_name in tissues_list:
        file_name = tissue_name + suffix
        file_path = os.path.join(path, file_name)
        print('Loading', file_name)

        with gzip.open(file_path) as file_object:
            header = {
                name: position
                for position, name in enumerate(
                    next(file_object).split()
                )
            }
            slope_pos = header['slope']
            variant_id_pos = header['variant_id']
            get_variant = itemgetter(variant_id_pos)
            get_slope = itemgetter(slope_pos)

            for line in tqdm(file_object):
                data = line.split()
                variant_id = get_variant(data)
                slope = get_slope(data)
                bdb[variant_id].add(tissue_name + ',' + slope)


"""For tests use:
from collections import namedtuple


def convert(dictionary):
    return namedtuple('GenericDict', dictionary.keys())(**dictionary)


v = convert({'chrom': 1, 'pos': 787151, 'ref': 'G', 'alt': 'A'})

bdb = ExpressionDatabase('expression_slope_in_tissues_by_mutation.db')
bdb.get_by_mutation(v)

print(bdb['1_693731_A_G_b37'])
"""
