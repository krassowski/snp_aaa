import json
from collections import defaultdict

from tqdm import trange

from biomart_data import BiomartData
from commands import command
from commands import SourceSubparser
from snp_parser import BiomartDataset
from cache import cacheable
from variant import Variant, BiomartVariant
from variant_sources import variants_getter


DEFAULT_BIOMART_URL = 'http://grch37.ensembl.org/biomart'


class VariantsData(BiomartData):

    def __init__(self, dataset=None, filters=None):

        if filters is None:
            filters = {}

        filters.update({
            'so_parent_name':
                [
                    'protein_altering_variant',   # inframe + frameshift
                    'synonymous_variant',
                    'missense_variant',
                    'frameshift_variant'
                ]
        })
        super(self.__class__, self).__init__(dataset, BiomartVariant, filters)


def gene_names_from_patacsdb_csv(limit_to=None):
    """
    Load gene names from given csv file. The default file was exported from
    sysbio.ibb.waw.pl/patacsdb database, from homo sapiens dataset.

    The gene names are expected to be located at second column in each row.
    The count of gene names to read can be constrained with use of an optional
    `limit_to` argument (default None denotes that all gene names
    from the file should be returned)
    """
    genes = set()
    count = 0
    with open('patacsdb_all.csv', 'r') as f:
        for line in f:
            if limit_to and count >= limit_to:
                break
            count += 1
            if line:
                data = [x.strip() for x in line.split(',')]
                genes.add(data[1])

    return list(genes)

biomart_args = SourceSubparser(
    'biomart',
    help='Arguments for biomart arguments fetching'
)

biomart_args.add_command(
    '--dataset',
    type=str,
    help='name of biomart dataset to be used, eg. hsapiens_snp',
    default='hsapiens_snp_som'
)

biomart_args.add_command(
    '--filters',
    type=json.loads,
    default={},
    help='Additional variants filters to be used when querying biomart'
)

biomart_args.add_command(
    '--genes_list',
    nargs='+',
    type=str,
    default=gene_names_from_patacsdb_csv(),
    help=(
        'list of human genes from which variants should be extracted for use in '
        'available analyses. By default all human genes from PATACSDB '
        'will be used. To use all human genes specify "all_human_genes".'
    )
)

biomart_args.add_command(
    '--biomart',
    type=str,
    help=(
        'URL of biomart to be used. Default: %s. '
        'To query GRCh38 use http://www.ensembl.org/biomart/ '
        'For mirrors of GRCh38 replace "www" with: "uswest", "useast" or "asia". '
        %
        DEFAULT_BIOMART_URL
    ),
    default=DEFAULT_BIOMART_URL
)


@biomart_args.command(
    '--step_size',
    '-s',
    type=int,
    default=25,
    help='When using list of genes to download variants, how many genes '
         'should be sent to biomart at once? Defaults to 25.'
)
def step_size(value, args):
    if value <= 0:
        print('Step size has to be positive!')


@variants_getter
def biomart(args):
    """Download variants from given --biomart, belonging to
    coding sequences of genes specified with --genes_list and
    filtered with --filters.
    """
    dataset = BiomartDataset(args.biomart, name=args.dataset)

    raw_variants_by_gene = download_variants.save(
        args.biomart,
        dataset,
        args.genes_list,
        step_size=args.step_size,
        filters=args.filters
    )
    return raw_variants_by_gene


@cacheable
def get_all_human_genes(biomart):
    """Retrieves list of human gene identifiers from Ensembl's biomart."""
    gene_names = set()

    print('Downloading list of all human genes in Ensembl...')

    genes_dataset = BiomartDataset(biomart, name='hsapiens_gene_ensembl')
    response = genes_dataset.search({'attributes': ['ensembl_gene_id']})

    for line in response.iter_lines():
        gene = line.decode('utf-8')
        gene_names.add(gene)

    print('Download of human genes list finished.')

    return gene_names


def find_difference(variant, known_variant):
    """Find differences between two Variants.

    Raises:
        ValueError if variants differ in anything which is not in __volatile_attributes__
    """
    variable_attributes = set(Variant.__volatile_attributes__)
    constant_attributes = set(Variant.__slots__) - variable_attributes

    if variant == known_variant:
        return True

    for attr in constant_attributes:
        if getattr(known_variant, attr) != getattr(variant, attr):
            raise ValueError(
                'Received variants with the same id: %s but different %s.'
                'Only: %s are allowed to differ in Biomart fetched data.'
                %
                (variant.refsnp_id, attr, ', '.join(variable_attributes))
            )


@cacheable
def download_variants(biomart, dataset, gene_names, step_size=50, filters={}):
    """Retrieve from Ensembl's biomart all variants that affect given genes.

    Variants that occur repeatedly in a given gene (i.e. were annotated for
    multiple transcripts) will be reported only once; list of affected
    transcripts will be preserved in 'affected_transcripts' property.

    The genes should be specified as the list of ensembl gene identifiers.


    Returns:
        dict of lists where variants are grouped by ensembl_gene_stable_id

        lists of variants are guaranteed to contain at least one variant
        each - if there is a key with id of given gene, it has at least one
        variant in list, fetched from Ensembl Biomart.
    """

    variants_by_gene = defaultdict(list)
    variants_count = 0
    transcripts_count = 0

    print('Downloading variants data from Ensembl\'s Biomart:')

    if gene_names == ['all_human_genes']:
        gene_names = get_all_human_genes.load_or_create(biomart)

    for start in trange(0, len(gene_names), step_size):

        query_filters = {'ensembl_gene': gene_names[start:start + step_size]}
        query_filters.update(filters)

        print('Downloading variants for genes:', gene_names[start:start + step_size])
        variants = VariantsData(filters=filters, dataset=dataset)
        print('Download completed. Parsing...')

        variants_in_gene = 0

        for biomart_variant in variants:
            assert biomart_variant.refsnp_id

            gene = biomart_variant.ensembl_gene_stable_id

            if gene not in gene_names:
                # print('Returned variant outside of query(!)')
                # print(gene_names)
                # print(biomart_variant)
                continue

            variants_with_the_same_id = [
                known_variant
                for known_variant in variants_by_gene[gene]
                if known_variant.refsnp_id == biomart_variant.refsnp_id
            ]

            assert len(variants_with_the_same_id) <= 1

            transcript = biomart_variant.extract_transcript()

            variant = Variant(
                affected_transcripts={transcript},
                **{
                    attr: getattr(biomart_variant, attr)
                    for attr in Variant.biomart_attributes
                }
            )
            transcripts_count += 1

            if variants_with_the_same_id:

                known_variant = variants_with_the_same_id[0]

                if variant != known_variant:
                    find_difference(variant, known_variant)

                known_variant.affected_transcripts.add(transcript)
            else:
                variants_by_gene[gene].append(variant)
                variants_in_gene += 1

        # this is a small trick to turn off unpickable iterator so it is
        # possible to save the variants object as a cache by pickling
        variants.iterator = None

        variants_count += variants_in_gene
        print('Parsed %s records.' % variants_in_gene)

    print('Downloaded %s transcripts and %s variants' % (transcripts_count, variants_count))
    return variants_by_gene


@command('--show_some_variant', action='store_true')
def show_some_variant(value, args):
    if not value:
        return

    variant = None
    dataset = BiomartDataset(args.biomart, name=args.dataset)

    genes_from_patacsdb = gene_names_from_patacsdb_csv()

    for gene in genes_from_patacsdb:

        variants_by_gene = download_variants.create(
            args.biomart,
            dataset,
            [gene]
        )

        if variants_by_gene:

            gene, variants = variants_by_gene.popitem()
            variant = variants[0]
            break

    if variant:
        print(variant)
    else:
        print('No variants found in given dataset.')

