from os.path import getsize, join

from numpy import log, max, sum
from pyLBL.arts_crossfit import CrossSection, download
import pytest


def test_artscrossfit(tmpdir, atmosphere, spectral_grid):
    download(tmpdir)
    formula = "CFC11"
    xsec = CrossSection(formula, join(tmpdir, "coefficients", f"{formula}.nc"))
    i = -1
    cross_section = xsec.absorption_coefficient(spectral_grid, atmosphere.t[i],
                                                atmosphere.p[i])
    assert log(max(cross_section)) == pytest.approx(-49.102292000921295)
    dv = spectral_grid[1] - spectral_grid[0]
    assert log(sum(cross_section)*dv) == pytest.approx(-46.04953693253788)


def test_download(tmpdir):
    dataset = {
        "C2F6.nc": 6374380,
        "CH3CCl3.nc": 3995343,
        "SF6.nc": 4158480,
        "CHCl3.nc": 4424122,
        "HFC125.nc": 4170602,
        "HCFC142b.nc": 7882125,
        "HCFC141b.nc": 15765917,
        "HFC143a.nc": 17705524,
        "HFC134a.nc": 28692364,
        "C8F18.nc": 2335741,
        "HFC227ea.nc": 4031126,
        "Halon1211.nc": 3928946,
        "CF4.nc": 4537000,
        "CFC113.nc": 4230484,
        "NF3.nc": 3928939,
        "HFC152a.nc": 6533367,
        "Halon1301.nc": 3988546,
        "CH2Cl2.nc": 3928942,
        "HFC245fa.nc": 155344,
        "CFC114.nc": 3680966,
        "CCl4.nc": 5110860,
        "cC4F8.nc": 5153911,
        "Halon2402.nc": 3962146,
        "CFC115.nc": 523723,
        "HCFC22.nc": 20721462,
        "HFC23.nc": 17535381,
        "CFC12.nc": 14987301,
        "C5F12.nc": 5857842,
        "SO2F2.nc": 3995341,
        "HFC365mfc.nc": 942146,
        "HFC32.nc": 5012929,
        "C4F10.nc": 3995341,
        "C3F8.nc": 3928940,
        "C6F14.nc": 2335741,
        "CFC11.nc": 5721655,
        "HFC4310mee.nc": 5957528,
        "HFC236fa.nc": 394544,
    }
    download(tmpdir)
    for key, value in dataset.items():
        assert value == getsize(join(tmpdir, "coefficients", key))
