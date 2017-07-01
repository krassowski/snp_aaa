from analyses import report, reporter
from snp_parser import all_poly_a_variants


@reporter
def summarize_poly_aaa_variants(variants_by_gene):

    variant_aaa_report = []

    for variant in all_poly_a_variants(variants_by_gene):

        for transcript in variant.affected_transcripts:

            if not transcript.poly_aaa:
                continue

            variant_aaa_report += [
                '\t'.join(map(str, [
                    variant.snp_id,
                    None,#variant.ensembl_gene_stable_id,
                    aaa_data.increased,
                    aaa_data.decreased,
                    aaa_data.change,
                    variant.ref,
                    alt,
                    transcript.ensembl_id,
                    transcript.cds_start,
                    transcript.cds_end,
                ]))
                for alt, aaa_data in transcript.poly_aaa.items()
            ]

    report(
        'poly aaa increase and decrease by variants',
        variant_aaa_report,
        [
            'snp_id', 'gene', 'poly_aaa_increase',
            'poly_aaa_decrease', 'poly_aaa_change', 'ref', 'alt', 'transcript'
        ]
    )
