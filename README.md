# pyLBL

This application aims to provide a simple, flexible framework for calculating line-by-line
spectra.

## Installation

#### HAPI2
This application depends on HAPI2.  As of now, the maintainers of HAPI2 have decided to
make that repository private.  *Thus only users who have been granted access to that
repository can install/run this application.*  We hope the maintainers will open access to
that repository to the public in the near future.  If you have 2-factor authentication
active for your account, cloning the HAPI2 repository via HTTPS requires a github
[personal access token](https://stackoverflow.com/questions/31305945/git-clone-from-github-over-https-with-two-factor-authentication).

#### Installing using pip
Make sure that you have the most recent version of `pip`, then run
the following command in the base directory of the repository:

```python
pip install .
```

## High-level API
The calculation of absorption coefficients requires a `Spectroscopy` object.  This object
will control the construction of a local spectral database and allow users to choose
which models are used to calculate the molecular lines, various molecular continua,
and absorption cross sections.  Currently, the supported models are as follows:

|component | models                                                                       |
|--------- | ---------------------------------------------------------------------------- |
|lines     | ["grtcode"](https://github.com/menzel-gfdl/pygrt/tree/grips-code)            |
|continua  | ["mt_ckd"](https://github.com/GRIPS-code/MT_CKD/tree/fortran-90-and-python)  |

For example, to create a `Spectroscopy` object using GRTcode to calculate the lines
and the MT-CKD continuum, use:

```python
from pyLBL import Spectroscopy

spectroscopy = Spectroscopy(hapi_config={"api_key": "<your HITRAN api key>",},
                            lines_backend="grtcode", continua_backend="mt_ckd")
```

#### Spectral database management
Spectral database managment is performed behind-the-schemes using HAPI2.  There are
several configurable options for HAPI2 that may be passed to the `Spectroscopy` constructor
in a dictionary, but the only one that is required is the `api_key`.   You must create
an account on the [HITRAN website](https://hitran.org) in order to get an api key.  A
more complete set of options includes:

```python
hapi_options = {
    "engine": "sqlite", # Type of database to create.
    "database": "local", # Name of the database (with a .db suffix added).
    "user": "root",
    "pass": null,
    "database_dir": ".", # Directory where the local database will be created.
    "debug": true,
    "echo": false,
    "display_fetch_url": false,
    "proxy": null,
    "host": "http://hitran.org", # Location of the remote database.
    "api_version": "v2",
    "tmpdir": "tmp", # Directory were temporary files will be created.
    "api_key": "<your HITRAN api key>", # HITRAN api key associated with your account.
    "info": "server_info.json"
}
```

#### User atmospheric inputs
Atmospheric inputs should be passed in an `xarray DataSet` object.  As an example,
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
             "layer": (["z",], [1,])
         },
    )
```

As shown in this example, the units of presure must be `Pa`, temperature must be `K`,
and mole fraction must be `mol mol-1`.  Users may define a dictionary specifying which
variables in the dataset should be read:

```python
mapping = {
    "play": "p", # name of pressure variable in dataset.
    "tlay": "t", # name of temperature variable in dataset.
    "mole_fraction: {
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

Spectral grid input should in wavenumber [cm-1] and be defined as a numpy array, for example:

```python
from numpy import arange
grid = arange(1., 5001., 0.1)
```

#### Absorption output
Absorption coefficients can be calculated using the input described above by running:

```python
absorption = spectroscopy.compute_absorption(self, atmosphere, grid, mapping=None)
```

The output is returned as an xarray `Dataset`.  An example of what this output looks
like in netCDF format is shown below.  For each molecule, the components of the full
spectra are split into separate indices along the `mechanism` dimension.  The full spectra
can be calculated by adding the different components together.

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
