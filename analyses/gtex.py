import sys

from analyses import report, reporter
from commands import AnalysisSubparser
from snp_parser import select_poly_a_related_variants
from snp_parser import all_poly_a_variants
from expression_database import ExpressionDatabase, ExpressedGenes, import_expressed_genes
from expression_database import import_expression_data


GTEX_DATABASE = 'expression_slope_in_tissues_by_mutation.db'
GTEX_GENES = 'expressed_genes.db'

gtex_args = AnalysisSubparser(
    'gtex',
    help='Management of GTEx cache'
)


@gtex_args.command('--reload', action='store_true')
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
def poly_aaa_vs_expression(variants_by_gene):

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

    aaa_variants_list = list(all_poly_a_variants(variants_by_gene))

    print(
        'Analysing %s poly_a related variants (out of %s total).'
        % (len(aaa_variants_list), len(variants_by_gene))
    )

    for variant in aaa_variants_list:

        for transcript in variant.affected_transcripts:

            if not transcript.poly_aaa:
                continue

            expression_data_by_alt = bdb.get_by_mutation(variant, transcript)

            transcript.expression = {}

            for alt, aaa_data in transcript.poly_aaa.iteritems():

                expression_data = expression_data_by_alt.get(alt, None)

                if not expression_data:
                    #print('No expression for', variant.refsnp_id, 'with alt:', alt)
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

        """
        gtex_report += [(
            sum('up' in v.expression.itervalues() for v in poly_a_related_variants),
            sum('down' in v.expression.itervalues() for v in poly_a_related_variants),
            sum(
                sum('up' == expr for expr in v.expression.itervalues())
                for v in poly_a_related_variants
            ),
            sum(
                sum('down' == expr for expr in v.expression.itervalues())
                for v in poly_a_related_variants
            ),
            sum(data.increased for v in poly_a_related_variants for data in v.poly_aaa.itervalues()),
            sum(data.decreased for v in poly_a_related_variants for data in v.poly_aaa.itervalues())
        )]
        """

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
