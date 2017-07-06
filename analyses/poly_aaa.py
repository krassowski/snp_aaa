from recordclass import recordclass
from analyses import report, reporter
from helpers import all_poly_a_variants


@reporter
def summarize_poly_aaa_variants(variants_by_gene):

    columns = [
        'snp_id', 'gene', 'poly_aaa_increase',
        'poly_aaa_decrease', 'poly_aaa_change',
        'chr', 'start', 'end', 'ref', 'alt',
        'transcript', 'cds_start', 'cds_end'
    ]

    Record = recordclass('RecordPolyA', columns)

    aaa_variants = []

    for variant in all_poly_a_variants(variants_by_gene):

        for transcript in variant.affected_transcripts:

            if not transcript.poly_aaa:
                continue

            for alt, aaa_data in transcript.poly_aaa.items():

                record = Record(
                    variant.snp_id,
                    None,#variant.ensembl_gene_stable_id    # TODO
                    aaa_data.increased,
                    aaa_data.decreased,
                    aaa_data.change,
                    variant.chr_name,
                    variant.chr_start,
                    variant.chr_end,
                    variant.ref,
                    alt,
                    transcript.ensembl_id,
                    transcript.cds_start,
                    transcript.cds_end
                )

                aaa_variants.append(record)

    report(
        'poly aaa increase and decrease by variants',
        aaa_variants,
        columns
    )

    print('Unique variants: %s' % len(aaa_variants))
    print('Variants identifiers: %s' % sum(v.snp_id.count(',') + 1 for v in aaa_variants))
