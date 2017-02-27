from analyses import report, reporter
from snp_parser import all_poly_a_variants


@reporter
def summarize_poly_aaa_variants(variants_by_gene):

    variant_aaa_report = []

    for variant in all_poly_a_variants(variants_by_gene):

        variant_aaa_report += [
            '\t'.join(map(str, [
                variant.refsnp_id,
                variant.ensembl_gene_stable_id,
                data.increased,
                data.decreased,
                data.change,
                alt
            ]))
            for alt, data in variant.poly_aaa.items()
        ]

    report(
        'poly aaa increase and decrease by variants',
        variant_aaa_report,
        [
            'snp_id', 'gene', 'poly_aaa_increase',
            'poly_aaa_decrease', 'poly_aaa_change', 'alt'
        ]
    )