from numpy import array_equal, log, max, sum
from pyLBL import Database, Spectroscopy
import pytest


def test_spectroscopy(atmosphere_dataset, spectral_grid, downloaded_database):
    database = Database(downloaded_database)
    spec = Spectroscopy(atmosphere_dataset, spectral_grid, database)
    molecules = spec.list_molecules()
    assert molecules[0] == "H2O"
    assert molecules[-1] == "HFC236fa"
    assert len(molecules) == 88


def test_absorption(single_layer_atmosphere, coarse_grid, downloaded_database):
    database = Database(downloaded_database)
    spec = Spectroscopy(single_layer_atmosphere, coarse_grid, database)
    beta = spec.compute_absorption(output_format="total")
    beta = beta.data_vars["absorption"]
    wavenumber = beta.coords["wavenumber"]
    assert max(beta.data) == pytest.approx(154.77712952851365)
    assert log(sum(beta.data)) == pytest.approx(7.212513759327571)
    assert beta.attrs["units"] == "m-1"
    assert array_equal(wavenumber.data, coarse_grid)
    assert wavenumber.attrs["units"] == "cm-1"


def test_spectroscopy_bad_lines_model(atmosphere_dataset, spectral_grid,
                                      downloaded_database):
    database = Database(downloaded_database)
    with pytest.raises(KeyError):
        _ = Spectroscopy(atmosphere_dataset, spectral_grid, database, lines_backend="foo")


def test_spectroscopy_bad_continua_model(atmosphere_dataset, spectral_grid,
                                         downloaded_database):
    database = Database(downloaded_database)
    with pytest.raises(KeyError):
        _ = Spectroscopy(atmosphere_dataset, spectral_grid, database,
                         continua_backend="foo")


def test_spectroscopy_bad_xsec_model(atmosphere_dataset, spectral_grid,
                                     downloaded_database):
    database = Database(downloaded_database)
    with pytest.raises(KeyError):
        _ = Spectroscopy(atmosphere_dataset, spectral_grid, database,
                         cross_sections_backend="foo")
