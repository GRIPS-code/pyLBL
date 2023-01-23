from glob import glob
from os.path import join
from setuptools import Extension, find_packages, setup


def c_gas_optics_lib():
    directory = "pyLBL/c_lib"
    src = ["{}/{}".format(directory, x) for x in 
           ["absorption.c", "spectra.c", "spectral_database.c", "voigt.c"]]
    return Extension("pyLBL.c_lib.libabsorption",
                     sources=src,
                     include_dirs=[directory,],
                     extra_compile_args=[],
                     extra_link_args=["-lsqlite3", "-lm"])


setup(
    name="pyLBL",
    version="0.0.0",
    author="R. Menzel, R. Pincus",
    author_email="",
    description="Line-by-line absorption calculators.",
    url="",
    python_requires=">=3.6",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: ",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        # Required depenedencies.
        "netCDF4",
        "numpy",
        "scipy",
        "sqlalchemy",
        "xarray",
        "mt_ckd @ git+http://github.com/GRIPS-code/MT_CKD@fortran-90-and-python",
        "arts_crossfit @ git+http://github.com/menzel-gfdl/arts-crossfit@make-package",

        # To build documentation.
        "Sphinx",
        "sphinxcontrib-apidoc",
        "sphinxcontrib-napoleon",
        "sphinx-autopackagesummary",

        # Other
        "pyarts",
    ],
    entry_points={
        "pyLBL" : ["Gas=pyLBL.c_lib.gas_optics:Gas",],
        "mt_ckd": ["CO2Continuum=mt_ckd.carbon_dioxide:CarbonDioxideContinuum",
                   "H2OForeignContinuum=mt_ckd.water_vapor:WaterVaporForeignContinuum",
                   "H2OSelfContinuum=mt_ckd.water_vapor:WaterVaporSelfContinuum",
                   "N2Continuum=mt_ckd.nitrogen:NitrogenContinuum",
                   "O2Continuum=mt_ckd.oxygen:OxygenContinuum",
                   "O3Continuum=mt_ckd.ozone:OzoneContinuum",
        ],
        "arts_crossfit": ["CrossSection=arts_crossfit.cross_section:CrossSection"],
        "arts": ["Gas=pyLBL.pyarts_frontend.frontend:PyArtsGas",],
    },
    ext_modules = [c_gas_optics_lib(),],
)
