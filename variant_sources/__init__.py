from helpers import decorator_maker

sources = {}

variants_source = decorator_maker(sources, 'variants loader')

from . import pkdb
from . import clinvar
from . import ensembl
