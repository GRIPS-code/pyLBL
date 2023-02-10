from setuptools import Extension, setup


def c_gas_optics_lib():
    """Defines c extension library."""
    directory = "pyLBL/c_lib"
    src = ["{}/{}".format(directory, x) for x in
           ["absorption.c", "spectra.c", "spectral_database.c", "voigt.c"]]
    return Extension("pyLBL.c_lib.libabsorption",
                     sources=src,
                     include_dirs=[directory,],
                     extra_compile_args=[],
                     extra_link_args=["-lsqlite3", "-lm"])


# Required dependencies.
install_requires = [
    "matplotlib",
    "netCDF4",
    "numpy",
    "scipy",
    "sqlalchemy",
    "xarray",
]


# Documentation dependencies.
doc_requires = [
    "sphinx",
    "sphinxcontrib-apidoc",
    "sphinxcontrib-napoleon",
    "sphinx-autopackagesummary",
    "sphinx_pangeo_theme",
]


# Optional dependencies.
extras_require = {
    "complete": install_requires,
    "docs": doc_requires,
    "arts": ["pyarts",]
}


# Entry points.
entry_points = {
   "pyLBL": ["Gas=pyLBL.c_lib.gas_optics:Gas",],
   "mt_ckd": [
       "CO2Continuum=pyLBL.mt_ckd.mt_ckd.carbon_dioxide:CarbonDioxideContinuum",
       "H2OForeignContinuum=pyLBL.mt_ckd.mt_ckd.water_vapor:WaterVaporForeignContinuum",
       "H2OSelfContinuum=pyLBL.mt_ckd.mt_ckd.water_vapor:WaterVaporSelfContinuum",
       "N2Continuum=pyLBL.mt_ckd.mt_ckd.nitrogen:NitrogenContinuum",
       "O2Continuum=pyLBL.mt_ckd.mt_ckd.oxygen:OxygenContinuum",
       "O3Continuum=pyLBL.mt_ckd.mt_ckd.ozone:OzoneContinuum",
   ],
   "arts_crossfit": ["CrossSection=pyLBL.arts_crossfit.arts_crossfit.cross_section:CrossSection"],
   "arts": ["Gas=pyLBL.pyarts_frontend.frontend:PyArtsGas",],
}


setup(
    name="pyLBL",
    version="0.0.1",
    description="Line-by-line absorption calculators.",
    url="https://github.com/GRIPS-code/pyLBL",
    author="pyLBL Developers",
    author_email="",
    license="LGPL-2.1",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: LGPL-2.1",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    packages=[
        "pyLBL",
        "pyLBL.c_lib",
        "pyLBL.webapi",
        "pyLBL.pyarts_frontend",
        "pyLBL.mt_ckd.mt_ckd",
        "pyLBL.arts_crossfit.arts_crossfit",
    ],
    install_requires=install_requires,
    extras_require=extras_require,
    entry_points=entry_points,
    ext_modules=[c_gas_optics_lib(), ],
    package_data={"": ["*.nc"], },
)
