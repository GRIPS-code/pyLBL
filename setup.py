from setuptools import find_packages, setup


setup(
    name="pyLBL",
    version="0.0.0",
    author="R. Menzel, R. Pincus",
    author_email="",
    description="Line-by-line absorption calculators.",
    url="",
    python_requires=">=3.5",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: ",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "pyarts",
        "pyrad @ git+http://github.com/menzel-gfdl/pylbl@add-continua",
        "pygrt @ git+http://github.com/menzel-gfdl/pygrt@main",
    ],
    entry_points={
        "arts": ["Gas=pyLBL.pyarts_frontend:PyArtsGas",],
        "grtcode" : ["Gas=pygrt.gas_optics:Gas",],
        "pyrad" : ["Gas=pyrad.optics.gas:Gas",
                   "CrossSection=pyrad.lbl.hitran.cross_sections:HitranCrossSection",
                   "CIA=pyrad.lbl.hitran.collision_induced_absorption:HitranCIA",
                   "H2OSelfContinuum=pyrad.lbl.continua.water_vapor:WaterVaporSelfContinuum",
                   "H2OForeignContinuum=pyrad.lbl.continua.water_vapor:WaterVaporForeignContinuum",
                   "O3Continuum=pyrad.lbl.continua.ozone:OzoneContinuum",
        ],
    },
)
