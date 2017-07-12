from collections import OrderedDict

from pathlib import Path

from helpers import decorator_maker


reports_dir = Path('reports')
reports_dir.mkdir(exist_ok=True)


def row_to_tsv(data):
    return '\t'.join(map(str, data))


def report(name, rows, column_names=()):
    """Generates list-based report quickly.

    File will be placed in 'reports' dir and name will be derived
    from 'name' of the report. The data should be an iterable
    collection of strings. Empty reports will not be created.
    """
    rows = list(rows)

    if not rows:
        print('Nothing to report for %s' % name)
        return

    if type(rows[0]) is not str:
        rows = [row_to_tsv(row) for row in rows]

    filename = name.replace(' ', '_') + '.tsv'
    path = reports_dir / filename

    with path.open('w') as f:
        if column_names:
            f.write('\t'.join(column_names) + '\n')
        f.write('\n'.join(rows))

    print('Created report "%s" with %s entries' % (name, len(rows)))

    return filename


analyses = OrderedDict()

reporter = decorator_maker(analyses, 'analysis')


from . import zcrb1
from . import cosmic
from . import spidex
from . import gtex
from . import poly_aaa
from . import gtex_vs_spidex
