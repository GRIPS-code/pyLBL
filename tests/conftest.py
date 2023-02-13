from collections import namedtuple
from ftplib import FTP
from os import environ
from os.path import isfile

from numpy import arange, asarray
import pytest
from xarray import Dataset


Atmos = namedtuple("Atmos", ["p", "t", "vmr"])


def variable(data, units, standard_name):
    """Create a tuple to aid in creating variables in an xarray Dataset.

    Args:
        data: Numpy array of data values.
        units: String physical units attribute.
        standard_name: String standard name attribute.

    Returns:
        A tuple containing the data needed to create a variable in an xarray Dataset.
    """
    return (["layer",], data, {"units": units, "standard_name": standard_name})


@pytest.fixture
def molecule_names():
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
    return names


@pytest.fixture
def spectral_grid():
    return arange(1., 3250., 0.1)


@pytest.fixture
def coarse_grid():
    return arange(1., 3000., 1.)


@pytest.fixture
def atmosphere(molecule_names):
    """Set conditions for a default test atmosphere.

    Returns:
        A Numpy array of pressure, temperature, and wavenumbers and a dictionary
        of molecule volume mixing ratios.
    """
    pressure = asarray([117., 1032., 11419., 98388.])  # [Pa].
    temperature = asarray([269.01, 227.74, 203.37, 288.99])  # [K].
    volume_mixing_ratio = {
        molecule_names["H2O"]: asarray([5.244536e-06, 4.763972e-06, 3.039952e-06,
                                        6.637074e-03]),
        molecule_names["CO2"]: asarray([0.00036, 0.00036, 0.00036, 0.00035999]),
        molecule_names["O3"]: asarray([2.936688e-06, 7.415223e-06, 2.609510e-07,
                                       6.859128e-08]),
        molecule_names["N2O"]: asarray([1.050928e-08, 1.319584e-07, 2.895416e-07,
                                        3.199949e-07]),
        molecule_names["CH4"]: asarray([2.947482e-07, 8.817705e-07, 1.588336e-06,
                                        1.700002e-06]),
        molecule_names["CO"]: asarray([3.621464e-08, 1.761450e-08, 3.315927e-08,
                                       1.482969e-07]),
        molecule_names["O2"]: asarray([0.209, 0.209, 0.2090003, 0.208996]),
        molecule_names["N2"]: asarray([0.78, 0.78, 0.78, 0.78]),
    }
    return Atmos(p=pressure, t=temperature, vmr=volume_mixing_ratio)


@pytest.fixture
def atmosphere_dataset(atmosphere):
    """Create an xarray Dataset for a test atmosphere.

    Args:
        pressure: Numpy array of pressure values [Pa].
        temperature: Numpy array of temperature values [K].
        volume_mixing_ratio: Dictionary of Numpy arrays of volume mixing ratios [mol mol-1].

    Returns:
        An xarray Dataset for a test atmosphere.
    """
    data_vars = {
       "pressure": variable(atmosphere.p, "Pa", "air_pressure"),
       "temperature": variable(atmosphere.t, "K", "air_temperature"),
    }
    for key, value in atmosphere.vmr.items():
        standard_name = f"mole_fraction_of_{key}_in_air"
        data_vars[key] = variable(value, "mol mol-1", standard_name)
    return Dataset(data_vars=data_vars)


@pytest.fixture
def single_layer_atmosphere(atmosphere):
    data_vars = {
       "pressure": variable(atmosphere.p[-1:], "Pa", "air_pressure"),
       "temperature": variable(atmosphere.t[-1:], "K", "air_temperature"),
    }
    for key, value in atmosphere.vmr.items():
        standard_name = f"mole_fraction_of_{key}_in_air"
        data_vars[key] = variable(value[-1:], "mol mol-1", standard_name)
    return Dataset(data_vars=data_vars)


@pytest.fixture
def downloaded_database():
    name = "pyLBL-2-7-23.db"
    if isfile(name):
        return name
    with FTP("ftp.gfdl.noaa.gov") as ftp:
        ftp.login()
        ftp.cwd(environ["FTP_DB_DIR"])
        ftp.retrbinary(f"RETR {name}", open(name, "wb").write)
    return name
