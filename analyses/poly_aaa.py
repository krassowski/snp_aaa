from collections import Counter

from recordclass import recordclass
from analyses import report, reporter
from helpers import all_poly_a_variants, IdMapper


class EnsemblTranscriptToProteinMapper(IdMapper):

    filename = 'ensembl_transcript_to_protein.gz'


@reporter
def summarize_poly_aaa_variants(variants):

    columns = [
        'snp_id', 'gene', 'poly_aaa_increase',
        'poly_aaa_decrease', 'poly_aaa_change',
        'chr', 'start', 'end', 'ref', 'alt',
        'transcript', 'cds_start', 'cds_end'
    ]
    #mapper = EnsemblTranscriptToProteinMapper()

    Record = recordclass('RecordPolyA', columns)

    aaa_records = []
    aaa_variants = set()
    up_variants = {}
    down_variants = {}
    cosmic_variants = []
    all_variants_ids = []
    variants_sources = Counter()
    transcripts = set()

    for variant in all_poly_a_variants(variants, preserve_sources=True):
        print(variant)
        if 'COSMIC' in variant.source:
            cosmic_variants.extend(variant.snp_id.split(','))
        all_variants_ids.extend(variant.snp_id.split(','))

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

                if aaa_data.increased:
                    up_variants[variant] = True
                if aaa_data.decreased:
                    down_variants[variant] = True
                transcripts.add(transcript.ensembl_id)

                aaa_records.append(record)

            aaa_variants.add(variant)

        for source in set(variant.source.split(',')):
            variants_sources[source] += 1

    report(
        'poly aaa increase and decrease by variants',
        aaa_records,
        columns
    )
    report(
        'poly aaa sources',
        variants_sources.items(),
        ['source', 'count']
    )
    report(
        'cosmic mutations id for querying',
        cosmic_variants
    )
    report(
        'all ids',
        all_variants_ids
    )

    print('Affected transcripts: %s' % len(transcripts))
    print('Down variants: %s' % len(down_variants))
    print('Up variants: %s' % len(up_variants))
    print('Unique variants: %s' % len(aaa_variants))
    print('Variants identifiers: %s' % sum(v.snp_id.count(',') + 1 for v in aaa_variants))
    print(variants_sources)
