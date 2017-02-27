import os
from collections import OrderedDict

from snp_parser import VERBOSITY_LEVEL


def report(name, data, column_names=()):
    """Generates list-based report quickly.

    File will be placed in 'reports' dir and name will be derived
    from 'name' of the report. The data should be an iterable
    collection of strings. Empty reports will not be created.
    """
    directory = 'reports'

    if not data:
        if VERBOSITY_LEVEL:
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


REPORTERS = OrderedDict()


def reporter(func):
    REPORTERS[func.__name__] = func
    return func

# Temporarily, TODO remove asterixes
from cosmic import *
from spidex import *
from gtex import *
from poly_aaa import *
