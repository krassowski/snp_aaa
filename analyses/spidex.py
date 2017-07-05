from __future__ import print_function
from collections import OrderedDict, defaultdict, Counter
from itertools import combinations
from operator import itemgetter

from tqdm import tqdm

from analyses import report, reporter
from scipy.stats import ks_2samp
from helpers import all_poly_a_variants
from settings import SPIDEX_LOCATION
from cache import cacheable
import tabix
from recordclass import recordclass
import numpy as np
from ggplot import ggplot, aes, geom_density, ggtitle, xlab, ylab, ggsave
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)


headers = [
    'chromosome', 'position', 'ref_allele', 'mut_allele',
    'dpsi_max_tissue', 'dpsi_zscore', 'gene', 'strand', 'transcript',
    'exon_number', 'location', 'cds_type', 'ss_dist', 'commonSNP_rs'
]

SpidexRecord = recordclass('SpidexRecord', headers)


def row_to_tsv(data):
    return '\t'.join(map(str, data))


def show_spanr_queries(to_test_online, step=40):
    """Prepare queries to be run on http://tools.genes.toronto.edu"""
    to_test_online = list(to_test_online)
    for start in range(0, len(to_test_online), step):

        query = [
            row_to_tsv(row)
            for row in to_test_online[start:start + step]
        ]
        query = '\n'.join(query)
        print('---')
        print(query)


def save_plot(plot, extension='svg', size=(19.2, 12), path='reports/', switch_to_next_figure=True):
    from matplotlib import axes
    import matplotlib as mpl
    mpl.rcParams['figure.figsize'] = '%s, %s' % size

    seaborns = [
        sns.axisgrid.JointGrid,
        sns.axisgrid.FacetGrid
    ]

    if type(plot) is ggplot:
        ggsave(filename=path + plot.title + '.' + extension, plot=plot, width=size[0], height=size[1])
    elif type(plot) in seaborns:
        plot.fig.set_size_inches(*size)
        plot.fig.savefig(path + plot.ax.title.get_text() + '.' + extension)
    elif type(plot) is axes.Subplot:
        plot.figure.set_size_inches(*size)
        plot.figure.savefig(path + plot.title.get_text() + '.' + extension)
    elif plot is plt:
        figure = plt.gcf()
        axes = plt.gca()
        figure.set_size_inches(*size)
        figure.savefig(path + axes.title.get_text() + '.' + extension)
    else:
        raise Exception('Unrecognized plot type: %s' % type(plot))
    if switch_to_next_figure:
        new_figure()


def prepare_data_frame(data_dict, melt=True):
    df = pd.DataFrame(OrderedDict(
        (key, pd.Series(value))
        for key, value in data_dict.items()
    ))
    if melt:
        df = pd.melt(df)
    return df


def spidex_get_variant(tb, variant):

    pos = [
        'chr' + variant.chr_name,
        variant.chr_start - 1,
        variant.chr_end
    ]

    try:
        records = [
            SpidexRecord(*record)
            for record in tb.query(*pos)
        ]
    except tabix.TabixError:
        print('TabixError:')
        print(pos, variant)
        return []

    return records


class StrandMismatch(Exception):
    pass


class Mismatch(Exception):
    pass


class Intronic(Exception):
    pass


class ToManyRecords(Exception):
    pass


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
            raise ToManyRecords(message)
        else:
            print(message, 'For variant:')
            print(variant)
            print('Those records have been received:')
            for r in relevant_records:
                print(r)

    record = relevant_records[0]

    if location and record.location == location:
        repr_data = (variant.snp_id, variant.chr_name, variant.chr_start)

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

    if variant.chr_strand != strand:
        repr_data = (variant.snp_id, variant.chr_strand, strand)
        message = (
            'Skipping record of id "%s": '
            'incorrect strand %s in variant vs %s in record' %
            repr_data
        )
        if strict:
            raise StrandMismatch(message, repr_data, variant)
        else:
            print(message)
            print(record.ref_allele)
            print(variant)
            return []

    if convert_to_strand(record.ref_allele, record.strand) != variant.ref:
        message = 'Reference mismatch for %s!' % variant.snp_id
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
    to_test_online = set()
    skipped_intronic = []
    skipped_strand_mismatch = []
    skipped_indels = []
    counter = Counter()

    all_columns = [
        'chr_name',
        'chr_start',
        'chr_end',
        'ref',
        'alt',
        'ensembl_gene_stable_id',
        'chr_strand',
        'ensembl_transcript_stable_id',
        'aaa_increased',
        'aaa_decreased',
        'aaa_change',
        'aaa_before',
        'aaa_after',
        'snp_id',
        'dpsi_max_tissue',
        'dpsi_zscore'
    ]

    Record = recordclass('SpidexRecord', all_columns)

    for variant in variants_list:
        counter['variants'] += 1

        records = spidex_get_variant(tb, variant)

        to_skip = []
        for alt in variant.alts:
            if variant.is_insertion(alt) or variant.is_deletion(alt):
                skipped_indels.append((variant, alt))
                to_skip.append(alt)

        for transcript in variant.affected_transcripts:

            if not transcript.poly_aaa:
                continue

            for alt, aaa_data in transcript.poly_aaa.items():

                if alt in to_skip:
                    continue

                variant_data = Record(
                    variant.chr_name,
                    variant.chr_start,
                    variant.chr_end,
                    variant.ref,
                    alt,
                    None,#variant.ensembl_gene_stable_id,
                    transcript.strand,
                    transcript.ensembl_id,
                    aaa_data.increased,
                    aaa_data.decreased,
                    aaa_data.change,
                    aaa_data.before,
                    aaa_data.after,
                    variant.snp_id,
                    None,
                    None
                )

                relevant_record = None

                try:
                    relevant_record = choose_record(
                        records,
                        variant,
                        alt,
                        convert_strands=True,
                        strict=True
                    )
                except StrandMismatch as e:
                    skipped_strand_mismatch.append(e.args[1])
                except Intronic as e:
                    skipped_intronic.append(e.args[1])
                except Mismatch:
                    counter['mismatch'] += 1
                except ToManyRecords:
                    counter['multiple_records'] += 1

                if not relevant_record:

                    if alt and len(variant.ref) == len(alt) == 1:
                        to_test_online.add(
                            (
                                variant.chr_name,
                                variant.chr_start,
                                variant.snp_id,
                                variant.ref,
                                alt
                            )
                        )

                else:
                    record = relevant_record

                    spidex_raw_report.append([variant.snp_id, alt, aaa_data, record])
                    spidex_raw_unique[
                        (
                            variant.chr_name,
                            variant.chr_start,
                            variant.chr_end,
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
                    ] = [variant.snp_id, alt, aaa_data, record]

                    variant_data.dpsi_max_tissue = record.dpsi_max_tissue
                    variant_data.dpsi_zscore = record.dpsi_zscore

                    spidex_report.append(variant_data)

    print(
        'Following mutations were nor found in SPIDEX'
        ' but may be found manually in SPANR'
    )
    show_spanr_queries(to_test_online)

    print('Skipped %s indels (SPIDEX does only 1-1 SNPs)' % len(skipped_indels))
    print('Analysed %s mutations.' % counter['variants'])
    print('Not indel and not SNV: %s' % counter['multiple_records'])
    print('Mismatches: %s' % counter['mismatch'])

    report(
        'spidex',
        spidex_report,
        all_columns
    )
    print('Accepted %s mutations' % len({record.snp_id for record in spidex_report}))
    print('Accepted %s unique mutations' % len({(record.chr_name, record.chr_start, record.chr_end, record.alt) for record in spidex_report}))

    report(
        'spidex_to_test_online',
        to_test_online,
        ['chr_name', 'chr_start', 'snp_id', 'ref', 'alt']
    )

    report(
        'spidex_skipped_intronic',
        skipped_intronic,
        ['snp_id', 'chr_name', 'chr_start']
    )

    report(
        'spidex_skipped_strand_mismatch',
        skipped_strand_mismatch,
        ['snp_id', 'chr_strand', 'SPIDEX_strand']
    )

    return spidex_raw_unique


def divide_variants_by_poly_aaa(raw_unique_report):

    spidex_raw_report = list(raw_unique_report.values())

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

    variants_groups = {
        'increase': variants_list(lambda aaa: aaa.increased),
        'decrease': variants_list(lambda aaa: aaa.decreased),
        'constant': variants_list(lambda aaa: aaa.change == 0),
        'all': variants_list(lambda aaa: True)
    }

    return variants_groups


def new_figure():
    figures = plt.get_fignums()
    if figures:
        plt.figure(max(figures) + 1)


def plot_aaa_vs_spidex(variants_groups):

    variants = variants_groups

    aaa_changes = sorted(set([x['change'] for x in variants['all']]))
    aaa_lengths = sorted(set([x['new_aaa_length'] for x in variants['all']]))

    def density_plot(by='dpsi_zscore', categorical=True):

        if categorical:
            data_dict = {
                'muts increasing AAA': np.array(
                    [x[by] for x in variants['increase']]
                ),
                'muts decreasing AAA': np.array(
                    [x[by] for x in variants['decrease']]
                ),
                'muts not changing AAA length': np.array(
                    [x[by] for x in variants['constant']]
                )
            }
        else:
            data_dict = OrderedDict(
                (change, np.array(
                    [
                        x[by]
                        for x in variants['all']
                        if x['change'] == change
                    ]
                ))
                for change in aaa_changes
                if len(
                    [
                        x[by]
                        for x in variants['all']
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
    save_plot(p)

    p = density_plot(by='max_dpsi')
    save_plot(p)

    p = density_plot(categorical=False)
    save_plot(p)

    data_dict = OrderedDict(
        (
            change,
            np.array([
                variant_data['dpsi_zscore']
                for variant_data in variants['all']
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
    save_plot(p)

    p = sns.lmplot(x='variable', y='value', data=df, x_jitter=0.25)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, observations visually jittered | change')
    p.ax.set_xlabel('AAA track length change resulting from given mutation [noised with visual jitter]')
    p.ax.set_ylabel('PSI z-score')
    save_plot(p)

    g = sns.boxplot(x='variable', y='value', data=df)
    g.axes.set_title('Boxplot: Poly AAA mutations and PSI z-score | change')
    g.set_xlabel('AAA track length change resulting from given mutation')
    g.set_ylabel('PSI z-score')
    save_plot(g)

    g = sns.violinplot(x='variable', y='value', data=df)
    g.axes.set_title('Violin: Poly AAA mutations and PSI z-score | change')
    g.set_xlabel('AAA track length change resulting from given mutation')
    g.set_ylabel('PSI z-score')
    save_plot(g)

    data_dict = OrderedDict(
        (
            length,
            np.array([
                variant_data['dpsi_zscore']
                for variant_data in variants['all']
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
    save_plot(p)

    p = sns.lmplot(x='variable', y='value', data=df, x_jitter=0.25)
    p.ax.set_title('Regression: Poly AAA mutations and PSI z-score, observations visually jittered | length')
    p.ax.set_xlabel('AAA length resulting from given mutation [noised with visual jitter]')
    p.ax.set_ylabel('PSI z-score')
    save_plot(p)

    g = sns.boxplot(x='variable', y='value', data=df)
    g.axes.set_title('Boxplot: Poly AAA mutations and PSI z-score | length')
    g.set_xlabel('AAA track length resulting from given mutation')
    g.set_ylabel('PSI z-score')
    save_plot(g)



@reporter
def poly_aaa_vs_spidex(variants_by_gene):
    """Analysis of poly A track changing mutations using data from SPIDEX."""
    aaa_variants_list = all_poly_a_variants(variants_by_gene)
    raw_report = spidex_from_list(aaa_variants_list)

    print('Unique points', len(raw_report))
    print('Plotting')
    variants_groups = divide_variants_by_poly_aaa(raw_report)
    plot_aaa_vs_spidex(variants_groups)
    print('ks test')
    spidex_aaa_ks_test(variants_groups, already_divided=True)


@reporter
def all_variants_vs_spidex(variants_by_gene):
    """The same as poly_aaa_vs_spidex but for all variants, not only poly(A) related."""
    raise NotImplementedError

    all_variants = []

    for gene_variants in variants_by_gene.values():
        all_variants.extend(gene_variants)

    # plot_aaa_vs_spidex(all_variants)


@cacheable
def get_all_zscore():
    """Get all dpsi_zscore from spidex database."""
    return _get_all_zscores()


def count_spidex():
    from multiprocess import count_lines
    return count_lines(SPIDEX_LOCATION)


def _get_all_zscores():
    zscores = []
    from multiprocess import fast_gzip_read

    print('Counting...')

    count = count_spidex()

    print('Loading...')

    with fast_gzip_read(SPIDEX_LOCATION) as f:
        header = next(f)
        get_dpsi_zscore = itemgetter(headers.index('dpsi_zscore'))
        for line in tqdm(f, total=count-1):
            try:
                data = line.rstrip('\n').split('\t')
                # record = SpidexRecord(*data)
                # zscores.append(record.dpsi_zscore)
                zscores.append(float(get_dpsi_zscore(data)))
            except Exception as e:
                print(e)
                continue

    return zscores


@reporter
def spidex_aaa_ks_test(variants_groups, already_divided=False):

    if not already_divided:
        aaa_variants_list = all_poly_a_variants(variants_groups)
        raw_report = spidex_from_list(aaa_variants_list)
        variants_groups = divide_variants_by_poly_aaa(raw_report)

    groups_zscores = {
        name: [point['dpsi_zscore'] for point in group]
        for name, group in variants_groups.items()
    }

    for group_1, group_2 in combinations(groups_zscores, 2):
        print('%s vs %s:' % (group_1, group_2))
        z_scores_1 = groups_zscores[group_1]
        z_scores_2 = groups_zscores[group_2]
        ks_result = ks_2samp(z_scores_1, z_scores_2)
        print(ks_result)

    groups_new_aaa_lengths = defaultdict(list)

    for name, group in variants_groups.items():
        for point in group:
            new_aaa_length = point['new_aaa_length']
            groups_new_aaa_lengths[new_aaa_length].append(point['dpsi_zscore'])

    group = None
    name = None

    ks_results = {}

    for new_aaa_length in sorted(groups_new_aaa_lengths):
        print(
            'All mutations causing poly_aaa to be <= %s vs all mutations causing poly_aaa to be > %s:'
            % (new_aaa_length, new_aaa_length)
        )
        z_scores_1 = [
            zscore
            for name, group in groups_new_aaa_lengths.items()
            for zscore in group
            if name <= new_aaa_length
        ]
        z_scores_2 = [
            zscore
            for name, group in groups_new_aaa_lengths.items()
            for zscore in group
            if name > new_aaa_length
        ]
        if not z_scores_2:
            print('No mutations causing poly_aaa to be > %s' % new_aaa_length)
            continue
        ks_result = ks_2samp(z_scores_1, z_scores_2)
        print(ks_result)
        ks_results[new_aaa_length] = - np.log(ks_result.pvalue)

    lengths = list(ks_results.keys())

    plt.hist(
        lengths,
        weights=list(ks_results.values()),
        bins=range(min(lengths), max(lengths)),
        rwidth=0.9
    )
    plt.xticks(lengths)

    plt.xlabel('Length of poly(A) track: $x$')
    plt.ylabel(r'$-\log($P-Value$)$')
    plt.title(
        'Ks-test for groups: '
        'mutations effecting in poly(A) length $\leq$ $x$ vs mutations effecting in poly(A) length > $x$'
    )
    plt.grid(True)
    save_plot(plt)


@reporter
def spidex_ks_test(variants_by_gene):

    tb = tabix.open(SPIDEX_LOCATION)

    full_spidex_zscore_dist = get_all_zscore.load_or_create()

    for gene, variant in variants_by_gene:

        variants_z_scores = []

        variants_z_scores.extend([
            float(record.dpsi_zscore)
            for record in spidex_get_variant(tb, variant)
         ])

        ks_result = ks_2samp(full_spidex_zscore_dist, variants_z_scores)

        print(gene)
        print(ks_result)
