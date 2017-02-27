from __future__ import print_function
from collections import OrderedDict
from snp_parser import report
from snp_parser import all_poly_a_variants
from snp_parser import SPIDEX_LOCATION
import tabix
from recordclass import recordclass
from ggplot import ggplot, aes, geom_density, ggtitle, xlab, ylab
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)

DRAW_PLOTS = True


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


def draw_plot(plot):
    if DRAW_PLOTS:
        from matplotlib import axes

        if type(plot) is ggplot:
            plot.draw().waitforbuttonpress()
        elif type(plot) is sns.axisgrid.FacetGrid:
            plot.fig.show()
        elif type(plot) is axes.Subplot:
            plot.figure.show()
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


def spidex_from_list(variants_list):

    headers = [
        'chromosome', 'position', 'ref_allele', 'mut_allele',
        'dpsi_max_tissue', 'dpsi_zscore', 'gene', 'strand', 'transcript',
        'exon_number', 'location', 'cds_type', 'ss_dist', 'commonSNP_rs'
    ]

    SpidexRecord = recordclass('SpidexRecord', headers)

    tb = tabix.open(SPIDEX_LOCATION)

    spidex_raw_report = []
    spidex_report = []
    to_test_online = []
    skipped_intronic = []
    skipped_strand_mismatch = []

    counter = 0

    for variant in variants_list:
        counter += 1

        pos = [
            'chr' + variant.chr_name,
            variant.chrom_start - 1,
            variant.chrom_end
        ]

        records = [
            SpidexRecord(*record)
            for record in tb.query(*pos)
        ]

        for alt, aaa_data in variant.poly_aaa.items():

            variant_data = [
                variant.chr_name,
                variant.chrom_start,
                variant.chrom_end,
                variant.ref,
                alt,
                variant.ensembl_gene_stable_id,
                variant.chrom_strand,
                variant.ensembl_transcript_stable_id,
                aaa_data.increased,
                aaa_data.decreased,
                aaa_data.change,
                aaa_data.before,
                aaa_data.after,
                variant.refsnp_id,
            ]

            relevant_records = [
                record
                for record in records
                if record.mut_allele == alt
            ]

            if not relevant_records:
                to_test_online.append(
                    [variant.chr_name, variant.chrom_start, variant.refsnp_id, variant.ref, alt or '-']
                )

            for record in relevant_records:

                if record.ref_allele != variant.ref:
                    print(
                        'Reference mismatch for %s!' %
                        variant.refsnp_id
                    )
                    continue
                if record.location == 'intronic':
                    repr_data = (variant.refsnp_id, variant.chr_name, variant.chrom_start)
                    print(
                        'Skipping intronic record for: %s from %s %s' %
                        repr_data
                    )
                    skipped_intronic.append(repr_data)
                    continue

                strand = '1' if record.strand == '+' else '0'

                if variant.chrom_strand != strand:
                    repr_data = (variant.refsnp_id, variant.chrom_strand, strand)
                    print(
                        'Skipping record for: %s - '
                        'incorrect strand %s in variant vs %s in record' %
                        repr_data
                    )
                    skipped_strand_mismatch.append(repr_data)
                    continue

                spidex_raw_report.append([variant, alt, aaa_data, record])

                record_data = variant_data
                record_data += [record.dpsi_max_tissue, record.dpsi_zscore]

                spidex_report.append(record_data)

    print(
        'Following mutations were nor found in SPIDEX'
        ' but may be found manually in SPANR'
    )
    show_spanr_queries(to_test_online)

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

    return spidex_raw_report


def plot_aaa_vs_spidex(spidex_raw_report):

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
                    [x[by] for x in variants_all if x['change'] == change]
                ))
                for change in aaa_changes
                if len([x[by] for x in variants_all if x['change'] == change]) > 1
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
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, estimator=mean')
    p.ax.set_xlabel('AAA track length change resulting from given mutation')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    p = sns.lmplot(x='variable', y='value', data=df, x_jitter=0.25)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, observations visually jittered')
    p.ax.set_xlabel('AAA track length change resulting from given mutation [noised with visual jitter]')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    plt.figure(max(plt.get_fignums()) + 1)
    g = sns.boxplot(x='variable', y='value', data=df)
    g.axes.set_title('Boxplot: Poly AAA mutations and PSI z-score')
    g.set_xlabel('AAA track length change resulting from given mutation')
    g.set_ylabel('PSI z-score')
    draw_plot(g)

    plt.figure(max(plt.get_fignums()) + 1)
    g = sns.violinplot(x='variable', y='value', data=df)
    g.axes.set_title('Violin: Poly AAA mutations and PSI z-score')
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
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, estimator=mean')
    p.ax.set_xlabel('AAA track length resulting from given mutation')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    p = sns.lmplot(x='variable', y='value', data=df, x_jitter=0.25)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, observations visually jittered')
    p.ax.set_xlabel('AAA length resulting from given mutation [noised with visual jitter]')
    p.ax.set_ylabel('PSI z-score')
    draw_plot(p)

    plt.figure(max(plt.get_fignums()) + 1)
    g = sns.boxplot(x='variable', y='value', data=df)
    g.axes.set_title('Boxplot: Poly AAA mutations and PSI z-score')
    g.set_xlabel('AAA track length resulting from given mutation')
    g.set_ylabel('PSI z-score')
    draw_plot(g)

    #from code import interact
    #interact(local=dict(globals(), **locals()))


def poly_aaa_vs_spidex(variants_by_gene):
    """Analysis of poly A track changing mutations using data from SPIDEX."""
    aaa_variants_list = all_poly_a_variants(variants_by_gene)
    raw_report = spidex_from_list(aaa_variants_list)
    plot_aaa_vs_spidex(raw_report)


def all_variants_vs_spidex(variants_by_gene):
    """The same as poly_aaa_vs_spidex but for all variants, not only poly(A) related."""
    all_variants = []
    for gene_variants in variants_by_gene.values():
        all_variants.extend(gene_variants)

    raw_report = spidex_from_list(all_variants)
    plot_aaa_vs_spidex(raw_report)
