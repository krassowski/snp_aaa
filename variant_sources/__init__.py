from analyses import decorator_maker

VARIANTS_GETTERS = {}

variants_getter = decorator_maker(VARIANTS_GETTERS, 'variants loader')

from . import pkdb
from . import clinvar
from . import ensembl
