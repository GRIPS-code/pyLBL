from ctypes import c_char_p, c_double, c_int, CDLL
from glob import glob
from pathlib import Path

from numpy import zeros
from numpy.ctypeslib import ndpointer


library = glob(str(Path(__file__).parent / "libabsorption*.so"))[0]
library = CDLL(library)


def check_return_code(value):
    if value != 0:
        raise ValueError("Error inside c functions.")
    return value


class Gas(object):
    """API for gas optics calculation."""
    def __init__(self, lines_database, formula):
        self.database = lines_database.path
        self.formula = formula

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid,
                               remove_pedestal=False, cut_off=25):
        v0 = int(round(grid[0]))
        vn = int(round(grid[-1]) + 1)
        n_per_v = int(round(1./(grid[1] - grid[0])))
        remove_pedestal = 1 if remove_pedestal else 0
        k = zeros((vn - v0)*n_per_v)
        library.absorption.argtypes = 3*[c_double,] + \
                                      3*[c_int,] + \
                                      [ndpointer(c_double, flags="C_CONTIGUOUS"),] + \
                                      2*[c_char_p,] + \
                                      2*[c_int,]
        library.absorption.restype = check_return_code
        library.absorption(
            c_double(pressure),
            c_double(temperature),
            c_double(volume_mixing_ratio),
            c_int(v0),
            c_int(vn),
            c_int(n_per_v),
            k,
            bytes(self.database, encoding="utf-8"),
            bytes(self.formula, encoding="utf-8"),
            c_int(cut_off),
            c_int(remove_pedestal),
        )
        return k
