from ftplib import FTP
from logging import getLogger
from os import mkdir
from os.path import join
from shutil import rmtree
from tarfile import TarFile
from uuid import uuid4

from netCDF4 import Dataset as ncDataset
from numpy import arange, asarray, zeros
from xarray import DataArray, Dataset

from pyLBL import Spectroscopy


info = getLogger(__name__).info


def _variable(data, units, standard_name):
    return (["layer",], data, {"units": units, "standard_name": standard_name})


class Circ(object):
    def __init__(self, case=1):
        """Creates CIRC xarray datasets.
        Args:
            case: Integer case number.
        """
        address = "ftp2.gfdl.noaa.gov"
        mb_to_pa = 100.
        tarfile = "/perm/GFDL_pubrelease/test_data/grtcode-data.tar.gz"
        names = {
            "H2O": "water_vapor",
            "CO2": "carbon_dioxide",
            "O3": "ozone",
            "N2O": "nitrous_oxide",
            "CO": "carbon_monoxide",
            "CH4": "methane",
            "O2": "oxygen",
            "N2": "nitrogen",
        }

        # Downloads CIRC data from GFDL's FTP site.
        self.tmp = str(uuid4())
        mkdir(self.tmp)
        archive = join(self.tmp, "archive.tar.gz")
        with open(archive, "wb") as tarball:
            with FTP(address) as ftp:
                ftp.login()
                info("Downloading {} from {}.".format(tarfile, address))
                ftp.retrbinary("RETR {}".format(tarfile), tarball.write)
        tarball = TarFile.open(name=archive, mode="r")
        path = join("grtcode-data", "circ", "circ-case{}.nc".format(case))
        tarball.extract(path, path=self.tmp)

        # Creates an xarray Dataset with the proper standard_name attributes.
        with ncDataset(join(self.tmp, path), "r") as dataset:
            temperature = dataset.variables["layer_temperature"][:]
            pressure = dataset.variables["layer_temperature"][:]*mb_to_pa
            vmr = {x: dataset.variables["{}_abundance".format(x)][:]
                   for x in ["H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2"]}
            vmr["N2" ] = asarray(pressure.size*[0.78])
            vars_ = {
                "play": _variable(pressure, "Pa", "air_pressure"),
                "tlay": _variable(temperature, "K", "air_temperature"),
            }
            for key, value in vmr.items():
                vars_.update({key: _variable(value, "mol mol-1",
                                             "mole_fraction_of_{}_in_air".format(names[key]))})
        self.dataset = Dataset(data_vars=vars_)

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        rmtree(self.tmp)


if __name__ == "__main__":
    with Circ(1) as circ:
        grid = arange(1., 3250., 1.)
        s = Spectroscopy(circ.dataset, grid,
                         hapi_config={"api_key": None})
        absorption_coefficient = s.compute_absorption(output_format="all")
        absorption_coefficient.to_netcdf("results.nc")
