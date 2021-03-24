from ftplib import FTP
from logging import getLogger
from os import mkdir
from os.path import join
from shutil import rmtree
from tarfile import TarFile
from tempfile import TemporaryDirectory
from uuid import uuid4

from numpy import arange, zeros
from xarray import DataArray, Dataset, open_dataset

from PyLBL import Spectroscopy


info = getLogger(__name__).info


class Circ(object):
    """Loads CIRC input data from a dataset.

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
        address = "ftp2.gfdl.noaa.gov"
        mb_to_pa = 100.
        tarfile = "/perm/GFDL_pubrelease/test_data/grtcode-data.tar.gz"

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
        data = {}
        with open_dataset(join(self.tmp, path)) as dataset:
            p = dataset["level_pressure"].copy()
            coord = {"levels": [i for i in range(p.data.size)]}
            data["pressure"] = p
            data["temperature"] = dataset["level_temperature"].copy()
            for molecule in ["CH4", "CO", "CO2", "H2O", "N2O", "O2", "O3", "N2"]:
                name = "vmr_{}".format(molecule); ab = "{}_abundance".format(molecule)
                data[name] = DataArray(zeros(p.data.size), dims=p.dims)
                if molecule == "N2":
                    data[name].values[:] = 0.78
                else:
                    data[name].values[1:] = dataset[ab][:]
                    data[name].values[0] = data[name].values[1]
        self.dataset = Dataset(data)
        self.dataset["pressure"].values[:] *= mb_to_pa

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        rmtree(self.tmp)

    @property
    def molecules(self):
        return [x.split("_")[-1] for x in self.dataset.data_vars.keys()
                if x.startswith("vmr_")]


if __name__ == "__main__":
    with Circ(1) as circ:
        grid = arange(1., 3250., 1.)
        spectroscopy = Spectroscopy(circ.molecules, "test.db")
        spectroscopy.load_spectral_inputs()
        absorption_coefficient = spectroscopy.compute_absorption(circ.dataset, grid)
        absorption_coefficient.to_netcdf("results.nc")
