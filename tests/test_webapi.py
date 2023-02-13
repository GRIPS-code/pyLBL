from collections import namedtuple
import pytest
from os import environ

from pyLBL import HitranWebApi, TipsWebApi
from pyLBL.webapi.hitran_api import NoIsotopologueError, NoTransitionsError
from pyLBL.webapi.tips_api import NoMoleculeError


Molecule = namedtuple("Molecule", ["id", "molecule_alias"])


def test_molecule_download():
    webapi = HitranWebApi(api_key=environ["HITRAN_API_KEY"])
    molecules = webapi.download_molecules()
    for i, id_, formula in zip([0, -1], [1, 1018], ["H2O", "Ar"]):
        assert molecules[i].id == int(id_)
        assert molecules[i].ordinary_formula == formula


def test_isotopologue_download():
    webapi = HitranWebApi(api_key=environ["HITRAN_API_KEY"])
    iso = webapi.download_isotopologues([Molecule(1, None),])
    for i, id_, molecule_id, name in zip([0, -1], [1, 129], [1, 1], ["H2", "D2"]):
        assert iso[i].id == id_
        assert iso[i].molecule_id == molecule_id
        assert iso[i].iso_name == f"{name}(16O)"


def test_transition_download():
    webapi = HitranWebApi(api_key=environ["HITRAN_API_KEY"])
    parameters = ["global_iso_id", "molec_id", "local_iso_id", "nu"]
    lines = webapi.download_transitions([Molecule(1, "H2O"),], 0, 3000, parameters)
    for i, id_, molecule_id, nu in zip([0, -1], [1, 1], [1, 1], [0.072049, 2999.90839]):
        assert lines[i].molec_id == molecule_id
        assert lines[i].local_iso_id == id_
        assert lines[i].nu == nu


def test_transition_download_no_iso():
    webapi = HitranWebApi(api_key=environ["HITRAN_API_KEY"])
    parameters = ["global_iso_id", "molec_id", "local_iso_id", "nu"]
    with pytest.raises(NoIsotopologueError):
        _ = webapi.download_transitions([], 0, 3000, parameters)


def test_transition_download_no_lines():
    webapi = HitranWebApi(api_key=environ["HITRAN_API_KEY"])
    parameters = ["global_iso_id", "molec_id", "local_iso_id", "nu"]
    with pytest.raises(NoTransitionsError):
        _ = webapi.download_transitions([Molecule(1, "H2O"),], 0, 1.e-12, parameters)


def test_tips_download():
    t, data = TipsWebApi().download("H2O")
    assert t[0] == 1.
    assert t[-1] == 6000.
    assert data[0, 0] == 1.
    assert data[-1, -1] == 1.1946e+07


def test_tips_download_no_molecule():
    with pytest.raises(NoMoleculeError):
        _, _ = TipsWebApi().download("XX")
