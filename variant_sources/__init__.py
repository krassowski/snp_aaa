VARIANTS_GETTERS = {}


def variants_getter(func):
    VARIANTS_GETTERS[func.__name__] = func
    return func

import biomart
import pkdb
import clinvar
