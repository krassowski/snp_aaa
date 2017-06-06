from analyses import decorator_maker

VARIANTS_GETTERS = {}

variants_getter = decorator_maker(VARIANTS_GETTERS, 'variants loader')

import pkdb
import clinvar
import ensembl
