from numpy import log, max
from pyLBL import Database, Gas
import pytest


def test_gas_optics(molecule_names, atmosphere, spectral_grid, downloaded_database):
    database = Database(downloaded_database)
    formula = "H2O"
    vmr_name = molecule_names[formula]
    layer = -1
    gas = Gas(database, formula)
    k = gas.absorption_coefficient(temperature=atmosphere.t[layer],
                                   pressure=atmosphere.p[layer],
                                   volume_mixing_ratio=atmosphere.vmr[vmr_name][layer],
                                   grid=spectral_grid)
    k = k[:spectral_grid.size]
    assert log(max(k)) == pytest.approx(-48.159224953962244)
    dv = spectral_grid[1] - spectral_grid[0]
    assert log(sum(k)*dv) == pytest.approx(-46.496121930910135)
