from numpy import sum
from pyLBL.mt_ckd.mt_ckd.carbon_dioxide import CarbonDioxideContinuum
from pyLBL.mt_ckd.mt_ckd.nitrogen import NitrogenContinuum
from pyLBL.mt_ckd.mt_ckd.oxygen import OxygenContinuum
from pyLBL.mt_ckd.mt_ckd.ozone import OzoneContinuum
from pyLBL.mt_ckd.mt_ckd.water_vapor import WaterVaporForeignContinuum
from pyLBL.mt_ckd.mt_ckd.water_vapor import WaterVaporSelfContinuum
import pytest


def create_vmr_dict(atmosphere, names, index):
    return {key: atmosphere.vmr[value][index] for key, value in names.items()}


def reference():
    values = {
        "CO2": [21.284607102488753, ],
        "H2OForeign": [131.87162317621952, ],
        "H2OSelf": [13.482864611247933, ],
        "N2": [0.7612890022253513, 0.5875825355004741, 0.00414557543788256, ],
        "O2": [0.24690308716508605, 0.11052072297118236, 0.03200556021322852,
               0.04514938962400228, 0.03897535512343981, 285.7607588975901,
               4419601.794329887, ],
        "O3": [0.0006562127133778276, 1.7334221226752753, 0.05197265302394795, ],
    }
    return values


def test_mt_ckd_co2(atmosphere, molecule_names):
    index = -1
    vmr = create_vmr_dict(atmosphere, molecule_names, index)
    values = reference()

    continua_dict = {
        "CO2": CarbonDioxideContinuum(),
        "H2OForeign": WaterVaporForeignContinuum(),
        "H2OSelf": WaterVaporSelfContinuum(),
        "N2":  NitrogenContinuum(),
        "O2": OxygenContinuum(),
        "O3": OzoneContinuum(),
    }

    for molecule, continua in continua_dict.items():
        for band, continuum in enumerate(continua.bands):
            spectra = continuum.spectra(atmosphere.t[index], atmosphere.p[index], vmr)
            assert values[molecule][band] == pytest.approx(sum(spectra))
