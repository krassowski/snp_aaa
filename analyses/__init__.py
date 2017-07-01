import os
from collections import OrderedDict

from helpers import decorator_maker


def report(name, data, column_names=()):
    """Generates list-based report quickly.

    File will be placed in 'reports' dir and name will be derived
    from 'name' of the report. The data should be an iterable
    collection of strings. Empty reports will not be created.
    """
    directory = 'reports'

    if not data:
        print('Nothing to report for %s' % name)
        return

    if not os.path.exists(directory):
        os.makedirs(directory)

    filename = name.replace(' ', '_') + '.txt'
    path = os.path.join(directory, filename)

    with open(path, 'w') as f:
        if column_names:
            f.write('\t'.join(column_names) + '\n')
        f.write('\n'.join(data))

    print('Created report "%s" with %s entries' % (name, len(data)))


analyses = OrderedDict()

reporter = decorator_maker(analyses, 'analysis')


from . import zcrb1
from . import cosmic
from . import spidex
from . import gtex
from . import poly_aaa
from . import gtex_vs_spidex
