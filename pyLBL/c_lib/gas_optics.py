"""Manages API or the molecular lines calcultion."""

from ctypes import c_char_p, c_double, c_int, CDLL
from glob import glob
from pathlib import Path

from numpy import zeros
from numpy.ctypeslib import ndpointer


library = glob(str(Path(__file__).parent / "libabsorption*.so"))[0]
library = CDLL(library)


def check_return_code(value):
    """Checks or errors occurring in the c routines.

    Args:
        value: Integer return code.

    Raises:
        ValueError if an error is encountered.
    """
    if value != 0:
        raise ValueError("Error inside c functions.")
    return value


class Gas(object):
    """API for gas optics calculation.

    Attributes:
        database: String path to the spectral sqlite3 database.
        formula: String chemical formula.
    """
    def __init__(self, lines_database, formula):
        """Initializes the object.

        Args:
            lines_database: Database object.
            formula: String chemical formula.
        """
        self.database = lines_database.path
        self.formula = formula

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid,
                               remove_pedestal=False, cut_off=25):
        """Calculates absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Numpy array defining the spectral grid [cm-1].
            remove_pedestal: Flag specifying if a pedestal should be subtracted.
            cut_off: Wavenumber cut-off distance [cm-1] from line centers.

        Returns:
            Numpy array of absorption coefficients [m2].
        """
        v0 = int(round(grid[0]))
        vn = int(round(grid[-1]) + 1)
        n_per_v = int(round(1./(grid[1] - grid[0])))
        remove_pedestal = 1 if remove_pedestal else 0
        k = zeros((vn - v0)*n_per_v)

        # Define argument types.
        library.absorption.argtypes = \
            3*[c_double,] + \
            3*[c_int,] + \
            [ndpointer(c_double, flags="C_CONTIGUOUS"),] + \
            2*[c_char_p,] + \
            2*[c_int,]

        # Set function to run on return.
        library.absorption.restype = check_return_code

        # Call the c function.
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
