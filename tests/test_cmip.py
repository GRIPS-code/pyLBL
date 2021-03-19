from glob import glob
from logging import getLogger
from os.path import join

from netCDF4 import Dataset
from numpy import arange, asarray, zeros
from xarray import DataArray, open_dataset

from PyLBL import Atmosphere, Spectroscopy


info = getLogger(__name__).info


class Cmip6(object):
    """Loads CIRC input data from a dataset.

    Attributes:
        pressure: Presssure [Pa] at layer interfaces.
        temperatuer: Temperature [K] at layer interfaces.
        layer_temperature: Temperature [K] at layer centers.
        surface_temperature: Surface temperature [K].
        vmr: Dictionary of volume-mixing ratios [mol mol-1] at layer interfaces.
    """
    def __init__(self, path):
        """Initializes object.

        Args:
            directory: Directory containing CMIP6 data files.
        """
        molecules = {"H2O": "hus", "CH4": "ch4", "CO": "co", "CO2": "co2", "N2O": "n2o",
                     "N2": "n2", "O2": "o2", "O3": "o3"}
        with open_dataset(path) as dataset:
            self.temperature = dataset["ta"].data[0,:,0,0]
            self.pressure = dataset["p"].data[0,:,0,0]
            self.vmr = {}
            for name, var in molecules.items():
                self.vmr[name] = dataset.variables[var].data[0,:,0,0]

if __name__ == "__main__":
    cmip = Cmip6("ta.nc")
    grid = arange(1., 3250., 1.)
    spectroscopy = Spectroscopy(cmip.vmr.keys(), "test.db")
    spectroscopy.load_spectral_inputs()
    print(spectroscopy.list_molecules())
    atmos = Atmosphere(cmip.temperature, cmip.pressure, cmip.vmr)
    absorption_coefficient = spectroscopy.compute_absorption(atmos, grid)
    with Dataset("results.nc", "w") as dataset:
        for name, units, data in zip(["pressure", "wavenumber"], ["Pa", "cm-1"],
                                     [atmos.pressure, grid]):
            dataset.createDimension(name, data.size)
            v = dataset.createVariable(name, "f8", (name,))
            v.setncattr("units", units)
            v[:] = data[:]
        for name, data in absorption_coefficient.items():
            v = dataset.createVariable("{}_absorption_coefficient".format(name), "f8",
                                       ("pressure", "wavenumber"))
            v.setncattr("units", "m-1")
            v[...] = data[...]
