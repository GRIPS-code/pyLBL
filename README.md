# pyLBL

This application aims to provide a simple, flexible framework for calculating line-by-line
spectra.

## Quickstart - calculating your first spectra
Let's jump right in with a very simple example.  After running through the installation
steps detailed below, spectra for a point in an atmosphere can be calculated by:

```python
from numpy import arange
from pyLBL import Database, Spectroscopy, WebApi
from xarray import Dataset

# Download the line paramters and store them in a local sqlite database.
# You should only have to create the database once.  Subsequent runs can
# then re-use the already created database.
webapi = WebApi("<your HITRAN api key>")
database = Database("<path to the database file>")
database.create(webapi)  # This step will take a long time.

# An xarray dataset describing the atmosphere to simulate is required.  This can be a
# cf-compliant netcdf file or explicitly defined (as below).
atmosphere = Dataset(
    data_vars={
        "play": (["z",], [98388.,], {"units": "Pa", "standard_name": "air_pressure"}),
        "tlay": (["z",], [288.99,], {"units": "K", "standard_name": "air_temperature"}),
        "xh2o": (["z",], [0.006637074,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_water_vapor_in_air"}),
        "xco2": (["z",], [0.0003599889,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_carbon_dioxide_in_air"}),
        "xo3": (["z",], [6.859128e-08,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_ozone_in_air"}),
        "xn2o": (["z",], [3.199949e-07,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_nitrous_oxide_in_air"}),
        "xco": (["z",], [1.482969e-07,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_carbon_monoxide_in_air"}),
        "xch4": (["z",], [1.700002e-06,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_methane_in_air"}),
        "xo2": (["z",], [0.208996,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_oxygen_in_air"}),
        "xn2": (["z",], [0.781,],
            {"units": "mol mol-1", "standard_name": "mole_fraction_of_nitrogen_in_air"}),
     },
)

# A wavenumber grid in units of [cm-1] is also required.
grid = arange(1., 5000., 0.1)

# A Spectroscopy object controls the absorption coefficient calculation.
s = Spectroscopy(atmosphere, grid, database)
spectra = s.compute_absorption()
spectra.to_netcdf("spectra.nc")
```

Here your HITRAN api key can be located by registering/logging into https://hitran.org and
looking at your profile page.  More information detailing each of the above steps are located
in the sections below.

## Installation

#### Installing using pip
Make sure that you have the most recent version of `pip`, then run
the following command in the base directory of the repository:

```python
pip install --upgrade pip # If you need to upgrade pip.
pip install .
```

This command should install the model and the following dependencies:

- matplotlib
- mt_ckd (https://github.com/GRIPS-code/MT_CKD/tree/fortran-90-and-python)
- netCDF4
- numpy
- scipy
- Sphinx
- sphinx-autopackagesummary
- sphinxcontrib-apidoc
- sphinxcontrib-napoleon
- SQLAlchemy
- xarray

## High-level API
This application aims to improve on existing line-by-line radiative transfer models
by separating the data management and calculation.  Data management is handled by
a `Database` object, which allows the user to construct a local database
of up-to-date line parameters, continuum coefficients, and cross sections
without having to explicitly interact with the ascii/data files typically
needed when using existing line-by-line models.  Absorption spectra calculation
is handled by a `Spectroscopy` object, which allows the user to specify which molecular
lines, continua, and cross section models they would like to use.  The details of
each of these objects are discussed further in the next sections.

#### Spectral database management
A `Database` object provides an interface to the current
[HITRAN](https://hitran.org) database of molecular line parameters.  To create a database
object of transitions for a specific set of molecules in a specific spectral range, run:

```python
from pyLBL import Database

# Make a connection to a database.  If the database already exists and you want to
# just to re-use it, this is the only step you need.
database = Database("<path to database>")

# If however you have not already populated the database, the data can be downloaded
# and inserted by running:
from pyLBL import WebApi
webapi = WebApi("<your HITRAN API key>")
database.create(webapi)  # Note that this step will take a long time.
```

You must create an account on the [HITRAN website](https://hitran.org) in order to get
an api key.  It is included as part of your profile on the webite.

#### Absorption calculation
A `Spectroscopy` object allow users to choose which models are used to calculate the
molecular lines, various molecular continua, and absorption cross sections.  Currently,
the supported models are as follows:

|component | models                                                                             |
|--------- | ---------------------------------------------------------------------------------- |
|lines     | ["pyLBL"](https://github.com/GRIPS-code/pyLBL/blob/new_db/pyLBL/spectral_lines.py) |
|continua  | ["mt_ckd"](https://github.com/GRIPS-code/MT_CKD/tree/fortran-90-and-python)        |

For example, to create a `Spectroscopy` object using the native pure python spectral lines model
and the MT-CKD continuum, use:

```python
from pyLBL import Spectroscopy

spectroscopy = Spectroscopy(atmosphere, grid, database, mapping=mapping,
                            lines_backend="pyLBL", continua_backend="mt_ckd")
```

Here the `database` argument is a `Database` object as described above.  The `atmosphere`,
`mapping`, and `grid` inputs are described in the following section.

#### User atmospheric inputs
Atmospheric inputs should be passed in as an xarray `Dataset` object.  As an example,
the surface layer of the first CIRC case can be described by:

```python
def variable(data, units, standard_name):
    return (["z",], data, {"units": units, "standard_name": standard_name})

def create_circ_xarray_dataset():
    from xarray import Dataset
    temperature = [288.99,] # [K].
    pressure = [98388.,] # [Pa].
    xh2o = [0.006637074,] # [mol mol-1].
    xco2 = [0.0003599889,] # [mol mol-1].
    xo3 = [6.859128e-08,] # [mol mol-1].
    xn2o = [3.199949e-07,] # [mol mol-1].
    xco = [1.482969e-07,] # [mol mol-1].
    xch4 = [1.700002e-06,] # [mol mol-1].
    xo2 = [0.208996,] # [mol mol-1].
    xn2 = [0.781,] # [mol mol-1].
    return Dataset(
        data_vars={
            "p": variable(pressure, "Pa", "air_pressure"),
            "t": variable(temperature, "K", "air_temperature"),
            "xh2o": variable(xh2o, "mol mol-1", "mole_fraction_of_water_vapor_in_air"),
            "xco2": variable(xco2, "mol mol-1", "mole_fraction_of_carbon_dioxide_in_air"),
            "xo3": variable(xo3, "mol mol-1", "mole_fraction_of_ozone_in_air"),
            "xn2o": variable(xn2o, "mol mol-1", "mole_fraction_of_nitrous_oxide_in_air"),
            "xco": variable(xco, "mol mol-1", "mole_fraction_of_carbon_monoxide_in_air"),
            "xch4": variable(xch4, "mol mol-1", "mole_fraction_of_methane_in_air"),
            "xo2": variable(xo2, "mol mol-1", "mole_fraction_of_oxygen_in_air"),
            "xn2": variable(xn2, "mol mol-1", "mole_fraction_of_nitrogen_in_air"),
         },
         coords={
             "layer": (["z",], [1,]),
         },
    )
```

As shown in this example, the units of presure must be Pa, temperature must be K,
and mole fraction must be mol mol<sup>-1</sup>.  Users may define a dictionary specifying which
variables in the dataset should be read:

```python
mapping = {
    "play": "p", # name of pressure variable in dataset.
    "tlay": "t", # name of temperature variable in dataset.
    "mole_fraction": {
        "H2O" : "xh2o", # name of water vapor mole fraction variable in dataset.
        "CO2" : "xco2", # name of carbon dioxided mole fraction variable in dataset.
        # et cetera
    },
}
```

If this dictionary is not provided, the application attempts to "discover" the variables
in the dataset using their CF `standard_name` attributes:

|variable                      | standard_name attribute         |
|----------------------------- | ------------------------------- |
|pressure                      | `"air_pressure"`                |
|temperature                   | `"air_temperature"`             |
|mole fraction of molecule xxx | `"mole_fraction_of_xxx_in_air"` |

For a full list of valid `standard_name` attributes, go [here](http://cfconventions.org/Data/cf-standard-names/77/build/cf-standard-name-table.html).

Spectral grid input should in wavenumber [cm<sup>-1</sup>] and be defined as a numpy
array, for example:

```python
from numpy import arange
grid = arange(1., 5001., 0.1)
```

#### Absorption output
Absorption coefficients can be calculated using the `Spectroscopy` object described
above by running:

```python
absorption = spectroscopy.compute_absorption(output_format="all")

# Optional: convert dataset to netcdf.
absorption.to_netcdf("<name of output file>")
```

The output is returned as an xarray `Dataset`.  The exact format of the output data
depends on the value of the `output_format` argument.  When set to `"all"` (which is
currently the default), the dataset will return the spectra split up by molecule
and mechansim (lines, continuum, cross_section). An example viewed in netCDF format
would look like this:

```
netcdf absorption {
dimensions:
        wavenumber = 49990 ;
        mechanism = 2 ;
        z = 1 ;
variables:
        double wavenumber(wavenumber) ;
                wavenumber:_FillValue = NaN ;
                wavenumber:units = "cm-1" ;
        string mechanism(mechanism) ;
        double H2O_absorption(z, mechanism, wavenumber) ;
                H2O_absorption:_FillValue = NaN ;
                H2O_absorption:units = "m-1" ;
        double CO2_absorption(z, mechanism, wavenumber) ;
                CO2_absorption:_FillValue = NaN ;
                CO2_absorption:units = "m-1" ;
        double O3_absorption(z, mechanism, wavenumber) ;
                O3_absorption:_FillValue = NaN ;
                O3_absorption:units = "m-1" ;
        double N2O_absorption(z, mechanism, wavenumber) ;
                N2O_absorption:_FillValue = NaN ;
                N2O_absorption:units = "m-1" ;
        double CO_absorption(z, mechanism, wavenumber) ;
                CO_absorption:_FillValue = NaN ;
                CO_absorption:units = "m-1" ;
        double CH4_absorption(z, mechanism, wavenumber) ;
                CH4_absorption:_FillValue = NaN ;
                CH4_absorption:units = "m-1" ;
        double O2_absorption(z, mechanism, wavenumber) ;
                O2_absorption:_FillValue = NaN ;
                O2_absorption:units = "m-1" ;
        double N2_absorption(z, mechanism, wavenumber) ;
                N2_absorption:_FillValue = NaN ;
                N2_absorption:units = "m-1" ;
data:

mechanism = "lines", "continuum" ;
}
```

If the `output_format` argument is instead set to `"gas"`, the spectra for
the different mechanims will be summed for each molecule, yielding output that looks
like this (in netCDF format):

```
netcdf absorption {
dimensions:
        wavenumber = 49990 ;
        z = 1 ;
variables:
        double wavenumber(wavenumber) ;
                wavenumber:_FillValue = NaN ;
                wavenumber:units = "cm-1" ;
        double H2O_absorption(z, wavenumber) ;
                H2O_absorption:_FillValue = NaN ;
                H2O_absorption:units = "m-1" ;
        double CO2_absorption(z, wavenumber) ;
                CO2_absorption:_FillValue = NaN ;
                CO2_absorption:units = "m-1" ;
        double O3_absorption(z, wavenumber) ;
                O3_absorption:_FillValue = NaN ;
                O3_absorption:units = "m-1" ;
        double N2O_absorption(z, wavenumber) ;
                N2O_absorption:_FillValue = NaN ;
                N2O_absorption:units = "m-1" ;
        double CO_absorption(z, wavenumber) ;
                CO_absorption:_FillValue = NaN ;
                CO_absorption:units = "m-1" ;
        double CH4_absorption(z, wavenumber) ;
                CH4_absorption:_FillValue = NaN ;
                CH4_absorption:units = "m-1" ;
        double O2_absorption(z, wavenumber) ;
                O2_absorption:_FillValue = NaN ;
                O2_absorption:units = "m-1" ;
        double N2_absorption(z, wavenumber) ;
                N2_absorption:_FillValue = NaN ;
                N2_absorption:units = "m-1" ;
}
```

Lastly, if the `output_format` argument is set to any other value, only the total
absorption spectra (summed over all molecules) will be returned.  In netCDF format, the
resulting dataset will appear like this:

```
netcdf absorption {
dimensions:
        wavenumber = 49990 ;
        z = 1 ;
variables:
        double wavenumber(wavenumber) ;
                wavenumber:_FillValue = NaN ;
                wavenumber:units = "cm-1" ;
        double absorption(z, wavenumber) ;
                absorption:_FillValue = NaN ;
                absorption:units = "m-1" ;
}
```
