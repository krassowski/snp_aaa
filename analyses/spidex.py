from __future__ import print_function
from collections import OrderedDict
import gzip

from tqdm import tqdm

from analyses import report, reporter
from snp_parser import jit
from snp_parser import all_poly_a_variants
from snp_parser import SPIDEX_LOCATION
from cache import cacheable
import tabix
from recordclass import recordclass
import numpy as np
from ggplot import ggplot, aes, geom_density, ggtitle, xlab, ylab, ggsave
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)

DRAW_PLOTS = True


headers = [
    'chromosome', 'position', 'ref_allele', 'mut_allele',
    'dpsi_max_tissue', 'dpsi_zscore', 'gene', 'strand', 'transcript',
    'exon_number', 'location', 'cds_type', 'ss_dist', 'commonSNP_rs'
]

SpidexRecord = recordclass('SpidexRecord', headers)


def row_to_tsv(data):
    return '\t'.join(map(str, data))


def show_spanr_queries(to_test_online, step=40, exclude_indels=True):
    """Prepare queries to be run on http://tools.genes.toronto.edu"""
    if exclude_indels:
        to_test_online = [
            [chrom, start, rs, ref, alt]
            for chrom, start, rs, ref, alt in to_test_online
            if alt != '-' and len(ref) == len(alt)
        ]
    for start in range(0, len(to_test_online), step):

        query = [
            row_to_tsv(row)
            for row in to_test_online[start:start + step]
        ]
        query = '\n'.join(query)
        print('---')
        print(query)


def draw_plot(plot, format='svg', size=(19.2, 12)):
    if DRAW_PLOTS:
        from matplotlib import axes
        import matplotlib as mpl
        mpl.rcParams['figure.figsize'] = '%s, %s' % size

        seaborns = [
            sns.axisgrid.JointGrid,
            sns.axisgrid.FacetGrid
        ]

        if type(plot) is ggplot:
            ggsave(filename=plot.title + '.' + format, plot=plot, width=size[0], height=size[1])
            #plot.draw().waitforbuttonpress()
        elif type(plot) in seaborns:
            #plot.fig.show()
            plot.fig.set_size_inches(*size)
            plot.fig.savefig(plot.ax.title.get_text() + '.' + format)
        elif type(plot) is axes.Subplot:
            #plot.figure.show()
            plot.figure.set_size_inches(*size)
            plot.figure.savefig(plot.title.get_text() + '.' + format)
        else:
            raise Exception('Unrecognized plot type: %s' % type(plot))


def prepare_data_frame(data_dict, melt=True):
    df = pd.DataFrame(OrderedDict(
        (key, pd.Series(value))
        for key, value in data_dict.iteritems()
    ))
    if melt:
        df = pd.melt(df)
    return df


def spidex_get_variant(tb, variant):

    pos = [
        'chr' + variant.chr_name,
        variant.chrom_start - 1,
        variant.chrom_end
    ]

    records = [
        SpidexRecord(*record)
        for record in tb.query(*pos)
    ]

    return records


class StrandMismatch(Exception):
    pass


class Mismatch(Exception):
    pass


class Intronic(Exception):
    pass


#@jit
def choose_record(records, variant, alt, location=None, convert_strands=False, strict=False):
    """
    location: if on of (None, 'intronic', 'exonic') is given,
        only variants from such location will be returned;
    strict: should exceptions be raised on anything
        suspicious or should we show just warnings instead?
    """

    relevant_records = []

    for record in records:
        if convert_to_strand(record.mut_allele, record.strand) in alt:
            relevant_records.append(record)

    if not relevant_records:
        return []

    if len(relevant_records) != 1:
        message = 'Too many records given!'
        if strict:
            raise ValueError(message)
        else:
            print(message, 'For variant:')
            print(variant)
            print('Those records have been received:')
            for r in relevant_records:
                print(r)

    record = relevant_records[0]

    if location and record.location == location:
        repr_data = (variant.refsnp_id, variant.chr_name, variant.chrom_start)

        message = (
            'Skipping intronic record for: %s from %s %s' %
            repr_data,
        )
        if strict:
            raise Intronic(message, repr_data)
        else:
            print(message)
            return []

    if convert_strands:
        strand = '1' if record.strand == '+' else '0'
    else:
        strand = record.strand

    if variant.chrom_strand != strand:
        repr_data = (variant.refsnp_id, variant.chrom_strand, strand)
        message = (
            'Skipping record for: %s - '
            'incorrect strand %s in variant vs %s in record' %
            repr_data
        )
        if strict:
            raise StrandMismatch(message, repr_data, variant)
        else:
            # TODO
            # print(message)
            # print(variant)
            return []

    if convert_to_strand(record.ref_allele, record.strand) != variant.ref:
        message = 'Reference mismatch for %s!' % variant.refsnp_id
        message += 'Variant ref: %s; record ref: %s, record strand: %s, record converted: %s' % (
            variant.ref,
            record.ref_allele,
            record.strand,
            convert_to_strand(record.ref_allele, record.strand)
        )
        if strict:
            raise Mismatch(
                message,
                convert_to_strand(record.ref_allele, record.strand),
                variant.ref,
                variant
            )
        else:
            print(message)
            print(variant)
            return []

    return record


complement = {
    'T': 'A',
    'A': 'T',
    'G': 'C',
    'C': 'G'
}


def convert_to_strand(sequence, strand):
    """If strand is -, return complementary sequence,
    else return the sequence unchanged"""
    if strand == '+':
        return sequence
    else:
        return ''.join([complement[nuc] for nuc in sequence])


def spidex_from_list(variants_list):

    tb = tabix.open(SPIDEX_LOCATION)

    spidex_raw_report = []
    spidex_raw_unique = {}
    spidex_report = []
    to_test_online = []
    skipped_intronic = []
    skipped_strand_mismatch = []
    skipped_indels = []

    counter = 0

    for variant in variants_list:
        counter += 1

        records = spidex_get_variant(tb, variant)

        to_skip = []
        for alt in variant.alts:
            if variant.is_insertion(alt) or variant.is_deletion(alt):
                skipped_indels.append((variant, alt))
                to_skip.append(alt)

        for transcript in variant.affected_transcripts:

            if not transcript.poly_aaa:
                continue

            for alt, aaa_data in transcript.poly_aaa.iteritems():

                if alt in to_skip:
                    continue

                variant_data = [
                    variant.chr_name,
                    variant.chrom_start,
                    variant.chrom_end,
                    variant.ref,
                    alt,
                    variant.ensembl_gene_stable_id,
                    transcript.strand,
                    transcript.ensembl_id,
                    aaa_data.increased,
                    aaa_data.decreased,
                    aaa_data.change,
                    aaa_data.before,
                    aaa_data.after,
                    variant.refsnp_id,
                ]

                relevant_record = None

                try:
                    #print('Records are:', records)
                    relevant_record = choose_record(
                        records,
                        variant,
                        alt,
                        convert_strands=True
                    )
                except StrandMismatch as e:
                    print(e.message)
                    skipped_strand_mismatch.append(e.args[1])
                except Intronic as e:
                    print(e.message)
                    skipped_intronic.append(e.args[1])
                except Mismatch as e:
                    print(e.message)

                if not relevant_record:
                    to_test_online.append(
                        [
                            variant.chr_name,
                            variant.chrom_start,
                            variant.refsnp_id,
                            variant.ref,
                            alt or '-'
                        ]
                    )

                else:
                    record = relevant_record

                    #print('Record', record)
                    spidex_raw_report.append([variant.refsnp_id, alt, aaa_data, record])
                    spidex_raw_unique[
                        (
                            variant.chr_name,
                            variant.chrom_start,
                            variant.chrom_end,
                            variant.ref,
                            alt,
                            transcript.strand,
                            aaa_data.increased,
                            aaa_data.decreased,
                            aaa_data.change,
                            aaa_data.before,
                            aaa_data.after,
                            record.dpsi_max_tissue,
                            record.dpsi_zscore
                        )
                    ] = [variant.refsnp_id, alt, aaa_data, record]

                    record_data = variant_data
                    #print('This record is of type ', record, ': >', record)
                    record_data += [record.dpsi_max_tissue, record.dpsi_zscore]

                    spidex_report.append(record_data)

    print(
        'Following mutations were nor found in SPIDEX'
        ' but may be found manually in SPANR'
    )
    show_spanr_queries(to_test_online)

    print('Skipped %s indels (SPIDEX does only 1-1 SNPs)' % len(skipped_indels))
    print('Analysed %s mutations.' % counter)

    report(
        'spidex',
        map(row_to_tsv, spidex_report),
        [
            'chr_name',
            'chrom_start',
            'chrom_end',
            'ref',
            'alt',
            'ensembl_gene_stable_id',
            'chrom_strand',
            'ensembl_transcript_stable_id',
            'aaa_increased',
            'aaa_decreased',
            'aaa_change',
            'aaa_before',
            'aaa_after',
            'refsnp_id',
            'dpsi_max_tissue',
            'record.dpsi_zscore'
        ]
    )
    report(
        'spidex_to_test_online',
        map(row_to_tsv, to_test_online),
        ['chr_name', 'chrom_start', 'refsnp_id', 'ref', 'alt']
    )

    report(
        'spidex_skipped_intronic',
        map(row_to_tsv, skipped_intronic),
        ['refsnp_id', 'chr_name', 'chrom_start']
    )

    report(
        'spidex_skipped_strand_mismatch',
        map(row_to_tsv, skipped_strand_mismatch),
        ['refsnp_id', 'chrom_strand', 'SPIDEX_strand']
    )

    return spidex_raw_unique


def plot_aaa_vs_spidex(spidex_raw_report, ignore_repetition=True):

    #if ignore_repetition:
    #    spidex_raw_report = list(set([(alt, aaa, rec) for v, alt, aaa, rec in spidex_raw_report]))
    #spidex_raw_report = list([[None, alt, aaa, SpidexRecord(*rec)] for alt, aaa, rec in spidex_raw_report])
    #print(spidex_raw_report)
    spidex_raw_report = list(spidex_raw_report.values())
    print('Unique points', len(spidex_raw_report))


    def variants_list(aaa_condition):
        return [
            {
                'new_aaa_length': aaa.after,
                'change': aaa.change,
                'max_dpsi': float(record.dpsi_max_tissue),
                'dpsi_zscore': float(record.dpsi_zscore)
            }
            for variant, alt, aaa, record in spidex_raw_report
            if aaa_condition(aaa)
        ]

    variants_increase = variants_list(lambda aaa: aaa.increased)
    variants_decrease = variants_list(lambda aaa: aaa.decreased)
    variants_constant = variants_list(lambda aaa: aaa.change == 0)
    variants_all = variants_list(lambda aaa: True)

    aaa_changes = sorted(set([x['change'] for x in variants_all]))
    aaa_lengths = sorted(set([x['new_aaa_length'] for x in variants_all]))

    def density_plot(by='dpsi_zscore', categorical=True):

        if categorical:
            data_dict = {
                'muts increasing AAA': np.array(
                    [x[by] for x in variants_increase]
                ),
                'muts decreasing AAA': np.array(
                    [x[by] for x in variants_decrease]
                ),
                'muts not changing AAA length': np.array(
                    [x[by] for x in variants_constant]
                )
            }
        else:
            data_dict = OrderedDict(
                (change, np.array(
                    [
                        x[by]
                        for x in variants_all
                        if x['change'] == change
                    ]
                ))
                for change in aaa_changes
                if len(
                    [
                        x[by]
                        for x in variants_all
                        if x['change'] == change
                    ]
                ) > 1
            )

        plot = (
            ggplot(
                aes(x='value', colour='variable', fill='variable'),
                data=prepare_data_frame(data_dict)
            ) +
            ggtitle('Impact of variants affecting poly AAA sequences on %s' % by) +
            xlab(by) +
            ylab('Kernel density estimate') +
            geom_density(alpha=0.6)
        )

        return plot

    p = density_plot()
    draw_plot(p)

    p = density_plot(by='max_dpsi')
    draw_plot(p)

    p = density_plot(categorical=False)
    draw_plot(p)

    data_dict = OrderedDict(
        (
            change,
            np.array([
                variant_data['dpsi_zscore']
                for variant_data in variants_all
                if variant_data['change'] == change
            ])
        )
        for change in aaa_changes
    )

    # seaborn
    # https://seaborn.pydata.org/tutorial/regression.html
    # http://seaborn.pydata.org/tutorial/categorical.html

    df = prepare_data_frame(data_dict)

    p = sns.lmplot(x='variable', y='value', data=df, x_estimator=np.mean)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, estimator=mean | change')
    p.ax.set_xlabel('AAA track length change resulting from given mutation')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    p = sns.lmplot(x='variable', y='value', data=df, x_jitter=0.25)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, observations visually jittered | change')
    p.ax.set_xlabel('AAA track length change resulting from given mutation [noised with visual jitter]')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    plt.figure(max(plt.get_fignums()) + 1)
    g = sns.boxplot(x='variable', y='value', data=df)
    g.axes.set_title('Boxplot: Poly AAA mutations and PSI z-score | change')
    g.set_xlabel('AAA track length change resulting from given mutation')
    g.set_ylabel('PSI z-score')
    draw_plot(g)

    plt.figure(max(plt.get_fignums()) + 1)
    g = sns.violinplot(x='variable', y='value', data=df)
    g.axes.set_title('Violin: Poly AAA mutations and PSI z-score | change')
    g.set_xlabel('AAA track length change resulting from given mutation')
    g.set_ylabel('PSI z-score')
    draw_plot(g)

    data_dict = OrderedDict(
        (
            length,
            np.array([
                variant_data['dpsi_zscore']
                for variant_data in variants_all
                if variant_data['new_aaa_length'] == length
            ])
        )
        for length in aaa_lengths
    )
    df = prepare_data_frame(data_dict)

    p = sns.lmplot(x='variable', y='value', data=df, x_estimator=np.mean)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, estimator=mean | length')
    p.ax.set_xlabel('AAA track length resulting from given mutation')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    p = sns.lmplot(x='variable', y='value', data=df, x_jitter=0.25)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, observations visually jittered | length')
    p.ax.set_xlabel('AAA length resulting from given mutation [noised with visual jitter]')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    plt.figure(max(plt.get_fignums()) + 1)
    g = sns.boxplot(x='variable', y='value', data=df)
    g.axes.set_title('Boxplot: Poly AAA mutations and PSI z-score | length')
    g.set_xlabel('AAA track length resulting from given mutation')
    g.set_ylabel('PSI z-score')
    draw_plot(g)

    #from code import interact
    #interact(local=dict(globals(), **locals()))


@reporter
def poly_aaa_vs_spidex(variants_by_gene):
    """Analysis of poly A track changing mutations using data from SPIDEX."""
    aaa_variants_list = all_poly_a_variants(variants_by_gene)
    raw_report = spidex_from_list(aaa_variants_list)
    plot_aaa_vs_spidex(raw_report)


@reporter
def all_variants_vs_spidex(variants_by_gene):
    """The same as poly_aaa_vs_spidex but for all variants, not only poly(A) related."""

    all_variants = []

    for gene_variants in variants_by_gene.itervalues():
        all_variants.extend(gene_variants)

    raise NotImplementedError
    # plot_aaa_vs_spidex(all_variants)


@cacheable
def get_all_zscore():
    """Get all dpsi_zscore from spidex database."""
    return _get_all_zscores()


@jit
def count_spidex():
    count = 0
    f = gzip.open(SPIDEX_LOCATION)
    for _ in f:
        count += 1
    f.close()
    return count


def _get_all_zscores():
    zscores = []

    print('Counting...')

    count = count_spidex()

    print('Loading...')

    f = gzip.open(SPIDEX_LOCATION)
    header = next(f)
    for line in tqdm(f, total=count-1):
        try:
            record = SpidexRecord(*line.rstrip('\n').split('\t'))
            zscores.append(record.dpsi_zscore)
        except Exception as e:
            print(e)
            continue
    f.close()

    return zscores


@reporter
def spidex_ks_test(variants_by_gene):
    """The same as poly_aaa_vs_spidex but for all variants, not only poly(A) related."""
    from scipy.stats import ks_2samp

    tb = tabix.open(SPIDEX_LOCATION)

    spidex_zscore_dist = get_all_zscore.load_or_create()

    for gene, variants in variants_by_gene:

        variants_z_scores = []

        for variant in variants:
            variants_z_scores.extend([
                record.dpsi_zscore
                for record in spidex_get_variant(tb, variant)
             ])

        ks_result = ks_2samp(spidex_zscore_dist, variants_z_scores)

        print(
            'Kolmogorov-Smirnov test for %s gene with %s variants analysed'
            %
            (gene, len(variants))
        )
        print(ks_result)
