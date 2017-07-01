from warnings import warn

try:
    from numba import jit
except ImportError:
    jit = lambda x: x
    warn('Install numba to speed up execution', UserWarning)
