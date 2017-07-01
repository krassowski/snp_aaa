from collections import defaultdict
from cna_by_transcript import CompleteCNA
from snp_parser import select_poly_a_related_variants
from snp_parser import COSMIC_VERSION
from analyses import report, reporter


def get_cosmic_genes_to_load(variants_by_gene):
    cosmic_gene_names = set()

    for gene, variants in variants_by_gene.items():
        cosmic_gene_names.update(
            [
                variant.gene
                for variant in variants
                ]
        )
    return cosmic_gene_names


def group_variants_for_cosmic(variants_by_gene):
    # COSMIC specific
    variants_by_gene_by_transcript = {}
    for gene, variants in variants_by_gene.items():
        by_transcript = defaultdict(list)

        for variant in variants:
            # The problem with the ensembl's biomart is that it returns records
            # from cosmic without information about the transcript, so we have
            # often a few identical records with only the snp_id different,
            # as for example: COSM3391893, COSM3391894
            # Fortunately the transcript id is encoded inside vcf_data retrieved
            # from biomart inside the gene identifier (if it is absent, then we
            # have a canonical transcript, at least it is the best guess), eg.:
            # ANKRD26_ENST00000376070 | ANKRD26_ENST00000436985 | ANKRD26
            gene_transcript_id = variant.gene
            transcript_id = gene_transcript_id

            by_transcript[transcript_id].append(variant)

        variants_by_gene_by_transcript[gene] = by_transcript
    return variants_by_gene_by_transcript


@reporter
def summarize_copy_number_expression(variants_by_gene):

    cosmic_genes_to_load = get_cosmic_genes_to_load(variants_by_gene)
    variants_by_gene_by_transcript = group_variants_for_cosmic(variants_by_gene)

    # TODO: expire cache on args.parse_variants, else always load.
    @cached(action='load')
    def cacheable_cna():
        return CompleteCNA(
            'cosmic/v' + COSMIC_VERSION + '/CosmicCompleteCNA.tsv',
            restrict_to=cosmic_genes_to_load
        )

    cna = cacheable_cna()

    variants_count = 0
    poly_a_related_variants_count = 0

    cnv_aaa_report = []
    import operator

    no_expression = set()

    for gene, variants_by_transcript in variants_by_gene_by_transcript.items():

        expression = [0, 0, 0]
        for gene_transcript_id in variants_by_transcript.keys():

            try:
                expr = cna.get_by_gene_and_transcript(gene_transcript_id)
                expression = map(operator.add, expression, expr)
            except KeyError:
                # no expression data for this transcript
                no_expression.add(gene_transcript_id)
                continue

        variants = []

        for transcript_variants in variants_by_transcript.values():
            variants.extend(transcript_variants)

        """
        # for this analysis consider only variants from cosmic
        variants = [
            variant
            for variant in variants
            if variant.source == 'COSMIC'
        ]
        """

        variants_count += len(variants)

        print(gene, 'with its', len(variants), 'variants analysed')

        # loss of poly_aaa, decrease in length (-1 per residue)
        decrease = 0
        # gain of poly_aaa: , increase in length (+1 per residue)
        increase = 0

        poly_a_related_variants = select_poly_a_related_variants(variants)
        poly_a_related_variants_count += len(poly_a_related_variants)

        # give scores for length increase
        for variant in poly_a_related_variants:
            if variant.poly_aaa.increased:
                increase += 1
            elif variant.poly_aaa.decreased:
                decrease += 1

        cnv_aaa_report += [
            (gene, expression[0], expression[2], increase, decrease)
        ]

    report('no expression data for some transcripts', no_expression)

    gene_count = len(variants_by_gene_by_transcript)

    print('analyzed genes:', gene_count)
    print('analyzed variants:', variants_count)
    print('poly_a related variants', poly_a_related_variants_count)

    report(
        'poly a and expression table',
        ['\t'.join(map(str, line)) for line in cnv_aaa_report],
        ['gene', 'cnv+', 'cnv-', 'aaa+', 'aaa-']
    )

