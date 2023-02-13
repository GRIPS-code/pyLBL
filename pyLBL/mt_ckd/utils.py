from os.path import dirname, join, realpath

from netCDF4 import Dataset
from numpy import asarray, copy, exp, interp, where, zeros


LOSCHMIDT = 2.6867775e19  # Loschmidt constant [cm-3]
P0 = 1013.25  # Reference pressure (1 atmosphere) [mb].
SECOND_RADIATION_CONSTANT = 1.4387752  # Second radiation constant [cm K].
T0 = 296.  # Reference temperature [K].
T273 = 273.15  # Reference temperature (0 celcius) [K].
m_to_cm = 100.  # [cm m-1].
Pa_to_mb = 0.01  # [mb Pa-1].


def air_number_density(pressure, temperature, volume_mixing_ratio):
    """Calculates the air number density.

    Args:
        pressure: Pressure [mb].
        temperature: Temperature [K].
        volume_mixing_ratio: Dictionary of volume mixing ratios [mol mol-1].

    Returns:
        Number density of air [cm-3].
    """
    return sum([dry_air_number_density(pressure, temperature, volume_mixing_ratio)*x
                for x in volume_mixing_ratio.values()])


def dry_air_number_density(pressure, temperature, volume_mixing_ratio):
    """Calculates the dry air number density.

    Args:
        pressure: Pressure [mb].
        temperature: Temperature [K].
        volume_mixing_ratio: Dictionary of volume mixing ratios [mol mol-1].

    Returns:
        Number density of dry air [cm-3].
    """
    return LOSCHMIDT*(pressure/P0)*(T273/temperature)*(1. - volume_mixing_ratio["H2O"])


def radiation_term(wavenumber, temperature):
    """Calculates the radiation term.

    Args:
        wavenumber: Array of wavenumber [cm-1].
        temperature: Temperature [K].

    Returns:
        The radiation term [cm-1].
    """
    t = temperature/SECOND_RADIATION_CONSTANT
    x = wavenumber[:]/t
    r = wavenumber[:]
    r = where(x <= 0.01, 0.5*x*wavenumber, r)
    return where(x <= 10., wavenumber*(1. - exp(-x))/(1. + exp(-x)), r)


def subgrid_bounds(grid, subgrid):
    """Calculates the starting and ending grid indices of a subgrid.

    Args:
        grid: A dictionary describing the main grid.
        subgrid: A dictionary describing the subgrid.

    Returns:
        The starting and ending grid indices of the subgrid.
    """
    if grid["resolution"] != subgrid["resolution"]:
        raise ValueError("grid and subgrid have different resolutions.")
    if grid["lower_bound"] > subgrid["lower_bound"] or \
       grid["upper_bound"] < subgrid["upper_bound"]:
        raise ValueError("subgrid not contained in grid.")
    lower = int((subgrid["lower_bound"] - grid["lower_bound"])/grid["resolution"])
    upper = int((subgrid["upper_bound"] - grid["lower_bound"])/grid["resolution"])
    return lower, upper


class Continuum(object):
    """Abstract class for gridded continuum coefficients."""
    def __init__(self, path):
        """Reads in the necessary data from an input dataset.

        Args:
            path: Path to the netcdf dataset.
        """
        raise NotImplementedError("You must override this class.")

    def spectra(self, temperature, pressure, vmr):
        """Calculates the spectral feature.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [mb].
            vmr: Dictionary of volume mixing ratios [mol mol-1].

        Return:
            An array of continuum extinction [cm-1].
        """
        raise NotImplementedError("You must override this method.")

    def grid(self):
        """Calculates the wavenumber grid [cm-1].

        Returns:
            A 1d numpy array containing the wavenumber grid [cm-1].
        """
        raise NotImplementedError("You must override this method.")


class Spectrum(object):
    """Helper class that reads data from a variable in the input dataset.

    Attributes:
        path: Path to the netcdf dataset.
        grid: Dictionary describing the wavenumber grid.
    """
    def __init__(self, path, name):
        """Reads the data from a variable in the input dataset.

        Args:
            path: Path to the netcdf dataset.
            name: Name of the variable in the dataset.
        """
        with Dataset(path, "r") as dataset:
            v = dataset.variables[name]
            self.data = copy(v[:])
            self.grid = {x: v.getncattr("wavenumber_{}".format(x)) for x in
                         ["lower_bound", "upper_bound", "resolution"]}
#           self.units = v.getncattr("units")

    def wavenumbers(self):
        """Calculates the wavenumber grid [cm-1] for the variable.

        Returns:
            A 1d numpy array containing the wavenumber grid [cm-1].
        """
        return asarray([self.grid["lower_bound"] + i*self.grid["resolution"]
                        for i in range(self.data.size)])


class BandedContinuum(object):
    """Contains all bands for a specific molecule's continuum.

    Attributes:
        bands: List of Continuum objects.
    """
    path = join(dirname(realpath(__file__)), "mt-ckd.nc")

    def __init__(self):
        """Reads in the necessary data from an input dataset."""
        raise NotImplementedError("You must override this method.")

    def spectra(self, temperature, pressure, vmr, grid):
        """Calculates the continum spectrum and interpolates to the input grid.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            vmr: Dictionary of volume mixing ratios [mol mol-1].
            grid: Array containing the spectral grid [cm-1].

        Return:
            An array of continuum extinction [m-1].
        """
        s = zeros(grid.size)
        for band in self.bands:
            s[:] += interp(grid, band.grid(),
                           band.spectra(temperature, pressure*Pa_to_mb, vmr),
                           left=0., right=0.)[:]*m_to_cm
        return s
