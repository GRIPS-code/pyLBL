from logging import getLogger
from os import mkdir
from os.path import join
from shutil import rmtree
from urllib.request import urlopen
from uuid import uuid4

from netCDF4 import Dataset as ncDataset
from numpy import arange, ones
from xarray import Dataset

from pyLBL import Database, Spectroscopy, WebApi


info = getLogger(__name__).info


def _variable(data, units, standard_name):
    return (["column", "layer"], data, {"units": units, "standard_name": standard_name})


class Rfmip(object):
    """Loads RFMIP input data from a dataset.

    Attributes:
        pressure: Presssure [Pa] at layer interfaces.
        temperatuer: Temperature [K] at layer interfaces.
        layer_temperature: Temperature [K] at layer centers.
        surface_temperature: Surface temperature [K].
        vmr: Dictionary of volume-mixing ratios [mol mol-1] at layer interfaces.
    """
    def __init__(self, case=1):
        """Initializes object.

        Args:
            case: Integer case number.
        """
        names = {"CH4": "methane_GM",
                 "CO": "carbon_monoxide_GM",
                 "CO2": "carbon_dioxide_GM",
                 "H2O": "water_vapor",
                 "N2": "nitrogen_GM",
                 "N2O": "nitrous_oxide_GM",
                 "O2": "oxygen_GM",
                 "O3": "ozone"}
        address = "/".join(["http://aims3.llnl.gov/thredds/fileServer/user_pub_work",
                            "input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/atmos/fx",
                            "multiple/none/v20190401",
                            "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"])
        self.tmp = str(uuid4()); path = "input.nc"
        mkdir(self.tmp)
        with open(join(self.tmp, path), "wb") as datafile:
            datafile.write(urlopen(address).read())

        with ncDataset(join(self.tmp, path), "r") as dataset:
            temperature = dataset.variables["temp_layer"][case - 1, ...]
            pressure = dataset.variables["pres_layer"]
            vmr = {}
            for key, value in names.items():
                units = float(dataset.variables[value].getncattr("units"))
                if value.endswith("_GM"):
                    vmr[key] = ones(pressure.shape)*dataset.variables[value][case - 1]*units
                else:
                    vmr[key] = dataset.variables[value][case - 1, ...]*units
            vars_ = {
                "play": _variable(pressure, "Pa", "air_pressure"),
                "tlay": _variable(temperature, "K", "air_temperature"),
            }
            for key, value in vmr.items():
                vars_.update({key: _variable(value, "mol mol-1",
                                             "mole_fraction_of_{}_in_air".format(names[key].strip("_GM")))})
            coords = {"column": ones(100), "layer": ones(60)}
            self.dataset = Dataset(data_vars=vars_, coords=coords)

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        rmtree(self.tmp)


if __name__ == "__main__":
    with Rfmip() as rfmip:
        webapi = WebApi(None)
        database = Database("foo.db")
        database.create(webapi, ["H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "N2"])
        grid = arange(1., 3250., 1.)
        spectroscopy = Spectroscopy(rfmip.dataset, grid, database)
        absorption_coefficient = spectroscopy.compute_absorption()
        absorption_coefficient.to_netcdf("results.nc")
