from ftplib import FTP
from logging import getLogger
from os.path import join
from tarfile import TarFile
from tempfile import TemporaryDirectory

from netCDF4 import Dataset
from numpy import arange, asarray, zeros

from PyLBL import Atmosphere, Spectroscopy


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
        with TemporaryDirectory() as tmp:
            archive = join(tmp, "archive.tar.gz")
            with open(archive, "wb") as tarball:
                with FTP(address) as ftp:
                    ftp.login()
                    info("Downloading {} from {}.".format(tarfile, address))
                    ftp.retrbinary("RETR {}".format(tarfile), tarball.write)
            tarball = TarFile.open(name=archive, mode="r")
            path = join("grtcode-data", "circ", "circ-case{}.nc".format(case))
            tarball.extract(path, path=tmp)
            with Dataset(join(tmp, path), "r") as dataset:
                self.pressure = dataset.variables["level_pressure"][:]*mb_to_pa
                self.temperature = dataset.variables["level_temperature"][:]
                self.layer_temperature = dataset.variables["layer_temperature"][:]
                self.surface_temperature = dataset.variables["surface_temperature"][:]
                self.vmr = {}
                for molecule in ["CH4", "CO", "CO2", "H2O", "N2O", "O2", "O3"]:
                    self.vmr[molecule] = zeros(self.pressure.size)
                    self.vmr[molecule][1:] = dataset.variables["{}_abundance".format(molecule)][:]
                    self.vmr[molecule][0] = self.vmr[molecule][1]
            self.vmr["N2"] = asarray(self.pressure.size*[0.78,])


if __name__ == "__main__":
    grid = arange(1., 3250., 1.)
    molecules = ["CH4", "CO", "CO2", "H2O", "N2", "N2O", "O2", "O3"]
    spectroscopy = Spectroscopy(molecules, "test.db")
    spectroscopy.load_spectral_inputs()
    print(spectroscopy.list_molecules())
    circ_data = Circ(1)
    atmos = Atmosphere(circ_data.temperature, circ_data.pressure, circ_data.vmr)
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
