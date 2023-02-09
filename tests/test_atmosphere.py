import pytest
from pyLBL.atmosphere import Atmosphere
from xarray import Dataset


def setup_atmosphere(names, dataset, mapping=None):
    atm = Atmosphere(dataset, mapping=mapping)
    for name in ["pressure", "temperature"]:
        assert getattr(atm, name).equals(dataset.data_vars[name])
    for key, value in names.items():
        assert atm.gases[key].equals(dataset.data_vars[value])


def test_atmosphere_without_mapping(molecule_names, atmosphere_dataset):
    setup_atmosphere(molecule_names, atmosphere_dataset)


def test_atmosphere_with_mapping(molecule_names, atmosphere_dataset):
    mapping = {
        "play": "pressure",
        "tlay": "temperature",
        "mole_fraction": molecule_names
    }
    setup_atmosphere(molecule_names, atmosphere_dataset, mapping)


def test_missing_standard_name(atmosphere_dataset):
    data_vars = dict(atmosphere_dataset.data_vars)
    del data_vars["pressure"]
    dataset = Dataset(data_vars=data_vars)
    with pytest.raises(ValueError):
        atm = Atmosphere(dataset)
