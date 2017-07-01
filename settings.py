GRCH_VERSION = 'GRCh37'
GRCH_SUBVERSION = '13'
ENSEMBL_VERSION = '88'
COSMIC_VERSION = '81'
DBSNP_VERSION = '150'
SPIDEX_LOCATION = 'spidex_public_noncommercial_v1.0/spidex_public_noncommercial_v1_0.tab.gz'
TRANSCRIPT_DB_PATH = 'ensembl/v' + ENSEMBL_VERSION + '/Homo_sapiens.' + GRCH_VERSION + '.cds.all.fa'
vcf_mutation_sources = {
    'COSMIC': {
        'is_alias': False,
        'path': 'cosmic/v' + COSMIC_VERSION + '/CosmicCodingMuts.vcf.gz.bgz',
        'given_as_positive_strand_only': True
    },
    'dbSNP': {
        'is_alias': False,
        'path': 'ncbi/dbsnp_' + DBSNP_VERSION + '-' + GRCH_VERSION.lower() + 'p' +
                GRCH_SUBVERSION + '/00-All.vcf.gz',
        'given_as_positive_strand_only': True
    },
    'ensembl': {
        'is_alias': False,
        'path': 'ensembl/v' + ENSEMBL_VERSION + '/Homo_sapiens.vcf.gz',
        'given_as_positive_strand_only': True
    },
    'ClinVar': {
        'is_alias': True,
        'aliased_vcf': 'dbSNP'
    },
    'ESP': {
        'is_alias': True,
        'aliased_vcf': 'ensembl'
    },
    'HGMD-PUBLIC': {
        'is_alias': True,
        'aliased_vcf': 'ensembl'
    },
    'PhenCode': {
        'is_alias': True,
        'aliased_vcf': 'ensembl'
    },
}
VERBOSITY_LEVEL = 0
