from os import environ
from os.path import normpath, sep

from pyLBL import Database, HitranWebApi


def test_database():
    molecules = ["H2O", "CO2", "CFC11"]
    database = Database(":memory:")
    webapi = HitranWebApi(api_key=environ["HITRAN_API_KEY"])
    database.create(webapi, molecules=molecules)

    db_molecules = database.molecules()
    for molecule in molecules:
        assert molecule in db_molecules

    formula, mass, transitions, _ = database.gas("H2O")
    assert formula == "H2O"
    assert mass[0] == 18.010565
    assert transitions[0].nu == 0.000134
    assert transitions[-1].nu == 41999.696489

    t, data = database.tips("CO2")
    assert t[0] == 1.
    assert t[-1] == 5000.
    assert data[0, 0] == 1.1722999811172485
    assert data[-1, -1] == 1018400.0

    path = database.arts_crossfit("CFC11")
    assert normpath(path).split(sep)[-3:] == [".cross-sections", "coefficients", "CFC11.nc"]
