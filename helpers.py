import datetime
import time
from collections import defaultdict

import seaborn as sns
from ggplot import ggplot, ggsave
from matplotlib import pyplot as plt

from jit import jit
from multiprocess import fast_gzip_read


def decorator_maker(storage, decorator_name):

    def decorator(func):
        name = func.__name__

        def informative_reporter(*args, **kwargs):
            print('Executing %s: %s...' % (decorator_name, name))
            with execution_time() as time:
                result = func(*args, **kwargs)
            print(
                'Execution of %s %s finished after %s'
                %
                (name, decorator_name, time.elapsed)
            )
            return result

        storage[name] = informative_reporter
        return informative_reporter

    return decorator


def select_poly_a_related_variants(variants):
    """Return list of variants occurring in a poly(A) track and variants which will create poly(A) track."""

    return [
        variant
        for variant in variants
        if any([
            data.has or data.will_have
            for affected_transcript in variant.affected_transcripts
            for data in affected_transcript.poly_aaa.values()
        ])
    ]


def group_variants(variants, preserve_sources=False):
    unique = {}
    for variant in variants:
        key = (variant.chr_name, variant.chr_start, variant.chr_end, variant.ref, variant.chr_strand, ''.join(sorted(variant.alts)))
        if key in unique:
            unique[key].snp_id += ',' + variant.snp_id
            unique[key].affected_transcripts = list(set(
                unique[key].affected_transcripts + variant.affected_transcripts
            ))
            if preserve_sources:
                unique[key].source += ',' + variant.source
        else:
            unique[key] = variant
    return unique.values()


def all_poly_a_variants(variants, merge_variants_with_multiple_id=True, preserve_sources=False):
    variants = variants.values()

    if merge_variants_with_multiple_id:
        variants = group_variants(variants, preserve_sources)

    poly_a_related_variants = select_poly_a_related_variants(variants)

    for variant in poly_a_related_variants:
        yield variant


class ExecutionTime:

    def __init__(self, start):
        self.start = start
        self.frozen = False

    @property
    def elapsed(self):
        if self.frozen:
            return self.frozen
        now = time.time()
        return datetime.timedelta(seconds=now - self.start)

    def freeze(self):
        self.frozen = self.elapsed


class execution_time:

    def __enter__(self):
        start = time.time()
        self.time = ExecutionTime(start)
        return self.time

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.time.freeze()


class IdMapper:

    filename = None

    def __init__(self, filename=None):
        data = defaultdict(list)
        if not filename:
            filename = self.filename
        if not filename:
            raise ValueError
        with fast_gzip_read(filename, processes=6) as f:
            header = next(f)
            assert header == '#hg19.knownToRefSeq.value	hg19.knownToEnsembl.value\n'
            for line in f:
                line = line.strip().split('\t')
                try:
                    ref_id, unknown_id = line
                except ValueError:
                    continue
                if unknown_id != 'n/a':
                    data[ref_id].append(unknown_id)
        self.data = data

    def map(self, ref_id):
        return self.data[ref_id]


@jit
def take_transcript_id_without_version(full_id):
    """Returns transcript id without version and everything which is after version separating comma.

    Example:
        Input: ESNT_0001.4
        Output: ESNT_0001

        Input: ESNT_0002.2.some_annotation
        Output: ESNT_0002

    """
    return full_id.split('.')[0]


def save_plot(plot, extension='svg', size=(10, 5), path='reports/', switch_to_next_figure=True, hide_title=False):
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
        title = axes.title.get_text()
        if hide_title:
            plt.title('')
        figure.set_size_inches(*size)
        figure.savefig(path + title + '.' + extension)
    else:
        raise Exception('Unrecognized plot type: %s' % type(plot))
    if switch_to_next_figure:
        new_figure()


def new_figure():
    figures = plt.get_fignums()
    if figures:
        plt.figure(max(figures) + 1)