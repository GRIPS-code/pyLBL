from ftplib import FTP
from logging import getLogger
from os import mkdir
from os.path import join
from shutil import rmtree
from tarfile import TarFile
from tempfile import TemporaryDirectory
from urllib.request import urlopen
from uuid import uuid4

from numpy import arange, zeros
from xarray import DataArray, Dataset, open_dataset

from PyLBL import Spectroscopy


info = getLogger(__name__).info


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
                 "O3": "ozone",
                 "pressure": "pres_layer",
                 "temperature": "temp_layer"}
        address = "/".join(["http://aims3.llnl.gov/thredds/fileServer/user_pub_work",
                            "input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/atmos/fx",
                            "multiple/none/v20190401",
                            "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"])

        self.tmp = str(uuid4()); path = "input.nc"
        mkdir(self.tmp)
        with open(join(self.tmp, path), "wb") as datafile:
            datafile.write(urlopen(address).read())
        data = {}
        with open_dataset(join(self.tmp, path)) as dataset:
            t = dataset[names["temperature"]]
            sizes = tuple([x for x in t.sizes.values()])
            for key, name in names.items():
                if name.endswith("_GM"):
                    varname = "vmr_{}".format(key)
                    data[varname] = DataArray(zeros(sizes), dims=t.dims)
                    for i in range(sizes[0]):
                        data[varname].values[i, ...] = dataset[name].data[i] * \
                                                       float(dataset[name].attrs["units"])
                elif key == "pressure":
                    data[key] = DataArray(zeros(sizes), dims=t.dims)
                    for i in range(sizes[0]):
                        data[key].values[i, ...] = dataset[name].data[...]
                elif key == "H2O" or key == "O3":
                    data["vmr_{}".format(key)] = dataset[name]
                else:
                    data[key] = dataset[name]
        self.dataset = Dataset(data)
        print(self.dataset)

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        rmtree(self.tmp)

    @property
    def molecules(self):
        return [x.split("_")[-1] for x in self.dataset.data_vars.keys()
                if x.startswith("vmr_")]


if __name__ == "__main__":
    with Rfmip() as rfmip:
        grid = arange(1., 3250., 1.)
        spectroscopy = Spectroscopy(rfmip.molecules, "test.db")
        spectroscopy.load_spectral_inputs()
        absorption_coefficient = spectroscopy.compute_absorption(rfmip.dataset, grid)
        absorption_coefficient.to_netcdf("results.nc")
