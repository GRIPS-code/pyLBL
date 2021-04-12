from glob import glob
from logging import getLogger

from numpy import arange
from xarray import Dataset, open_dataset

from pyLBL import Spectroscopy


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
        data = {}
        with open_dataset(path) as dataset:
            data["temperature"] = dataset["ta"]
            data["pressure"] = dataset["p"]
            for name, var in molecules.items():
                data["vmr_{}".format(name)] = dataset.variables[var]
        self.dataset = Dataset(data)

    @property
    def molecules(self):
        return [x.split("_")[-1] for x in self.dataset.data_vars.keys()
                if x.startswith("vmr_")]

if __name__ == "__main__":
    cmip = Cmip6("tests/ta.nc")
    grid = arange(1., 3250., 1.)
    spectroscopy = Spectroscopy(cmip.molecules, "test.db")
    spectroscopy.load_spectral_inputs()
    absorption_coefficient = spectroscopy.compute_absorption(cmip.dataset, grid)
    absorption_coefficient.to_netcdf("results.nc")
