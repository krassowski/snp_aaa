import sys

from analyses import report, reporter
from commands import command
from snp_parser import select_poly_a_related_variants
from snp_parser import all_poly_a_variants
from expression_database import ExpressionDatabase, ExpressedGenes, import_expressed_genes
from expression_database import import_expression_data


GTEX_DATABASE = 'expression_slope_in_tissues_by_mutation.db'
GTEX_GENES = 'expressed_genes.db'


@command('--reload_gtex', action='store_true')
def reload_gtex(value, args):
    if not value:
        return

    print('Reloading GTEx expression data:')

    bdb = ExpressedGenes(GTEX_GENES)
    bdb.reset()

    import_expressed_genes(
        bdb,
        path='GTEx_Analysis_v6p_eQTL',
        suffix='_Analysis.v6p.egenes.txt.gz'
    )

    bdb = ExpressionDatabase(GTEX_DATABASE)
    bdb.reset()

    import_expression_data(
        bdb,
        path='GTEx_Analysis_v6p_eQTL',
        suffix='_Analysis.v6p.signif_snpgene_pairs.txt.gz'
    )


@reporter
def gtex_over_api(variants_by_gene):
    import requests

    for variant in all_poly_a_variants(variants_by_gene):

        server = 'http://rest.ensembl.org'
        ext = '/eqtl/variant_name/homo_sapiens/%s?statistic=p-value;content-type=application/json' % variant.refsnp_id

        try:
            r = requests.get(
                server + ext,
                headers={
                    'content-type': 'application/json'
                }
            )

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            decoded = r.json()

            if not decoded['error']:
                print('Found sth!')
                print(repr(decoded))

        except Exception:
            pass


@reporter
def poly_aaa_vs_expression(variants_by_gene, include_all=False):

    bdb = ExpressionDatabase(GTEX_DATABASE)

    def is_length_difference_big(l1, l2):
        """Is the first list much longer than the second?"""
        len1 = len(l1)
        len2 = len(l2)
        assert len1 > len2

        if len2 == 0 or len1 / len2 > 10:
            return True

    gtex_report = []
    gtex_report_by_genes = []

    for gene, variants in variants_by_gene.iteritems():

        poly_a_related_variants = select_poly_a_related_variants(variants)

        print(
            'Analysing %s poly_a related variants (total: %s) from %s gene.'
            % (len(poly_a_related_variants), len(variants), gene)
        )

        for variant in poly_a_related_variants:

            variant.expression = {}
            expression_data_by_alt = bdb.get_by_mutation(variant)

            for alt, expression_data in expression_data_by_alt.items():

                if not expression_data:
                    print('No expression for', variant.refsnp_id)
                    if not include_all:
                        continue
                else:
                    print('Expression data for', variant.refsnp_id, 'found:', expression_data)

                expression_up = []
                expression_down = []

                for tissue, slope in expression_data:
                    if slope > 0:
                        expression_up += [tissue]
                    elif slope < 0:
                        expression_down += [tissue]

                # is this rather up?
                if len(expression_up) > len(expression_down):
                    # is this certainly up?
                    if is_length_difference_big(expression_up, expression_down):
                        expression_trend = 'up'
                    else:
                        expression_trend = 'rather_up'
                # is this rather down?
                elif len(expression_down) > len(expression_up):
                    # is this certainly down?
                    if is_length_difference_big(expression_down, expression_up):
                        expression_trend = 'down'
                    else:
                        expression_trend = 'rather_down'
                # is unchanged?
                else:
                    expression_trend = 'constant'

                expression_up_in_x_cases = len(expression_up)
                expression_down_in_x_cases = len(expression_down)

                data = variant.poly_aaa[alt]
                variant.expression[alt] = expression_trend

                gtex_report += [(
                    variant.refsnp_id,
                    expression_up_in_x_cases,
                    expression_down_in_x_cases,
                    expression_trend,
                    data.increased,
                    data.decreased,
                    data.change,
                    alt
                )]

        gtex_report_by_genes += [(
            gene,
            sum('up' in v.expression.values() for v in poly_a_related_variants),
            sum('down' in v.expression.values() for v in poly_a_related_variants),
            sum(
                sum('up' == expr for expr in v.expression.values())
                for v in poly_a_related_variants
            ),
            sum(
                sum('down' == expr for expr in v.expression.values())
                for v in poly_a_related_variants
            ),
            sum(data.increased for v in poly_a_related_variants for data in v.poly_aaa.values()),
            sum(data.decreased for v in poly_a_related_variants for data in v.poly_aaa.values())
        )]

    report(
        'Expression table for variants (based on data from GTEx)',
        ['\t'.join(map(str, line)) for line in gtex_report],
        [
            'variant', 'expression+', 'expression-', 'trend',
            'aaa+', 'aaa-', 'aaa change'
        ]
    )

    report(
        'Expression table for genes (based on data from GTEx)',
        ['\t'.join(map(str, line)) for line in gtex_report_by_genes],
        # note: alleles is not the same as variants
        [
            'gene', 'alleles with expression+', 'alleles with expression-',
            'variants with expression+', 'variants with expression-', '#aaa+', '#aaa-'
        ]
    )

    print('Done')
