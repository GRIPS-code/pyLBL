Atmospheric Inputs
------------------

This package aims to integrate easily into the pangeo suite of applications.  In view of
that, atmospheric inputs are required to be passed in as xarray :code:`Dataset` objects.
Along with the :code:`Dataset`, a user can pass in a dictionary with maps fixed names to
the variable names in the file.  The required keys are:

.. code-block:: python

  user_dictionary = {
      "play": "<name of the pressure variable in the dataset>",
      "tlay": "<name of the temperature variable in the dataset>",
      "mole_fraction": {
          "H2O": "<name of the water vapor variable in the dataset>",
          "CO2": "<name of the carbon dioxide variable in the dataset>",
          "O3": "<name of the ozone variable in the dataset>",
          "CH4": "<name of the methane variable in the dataset>",
          "N2O": "<name of the nitrous oxide variable in the dataset>",
          "CO": "<name of the carbon monoxide variable in the dataset>",
          # Other molecules of interest if desired.
      }
  }

If the user does not pass in the dictionary, the application attempts to "discover" the
correct variables in the :code:`Dataset` by examining the variables' CF
:code:`standard_name` attributes:

============================= ============================= ==============
Variable                      standard_name Attribute       Expected Units
============================= ============================= ==============
pressure                      "air_pressure"                "Pa"
temperature                   "air_temperature"             "K"
mole fraction of molecule xxx "mole_fraction_of_xxx_in_air" "mol mol-1"
============================= ============================= ==============

In either case, if the variables in the :code:`Dataset` do not have the expected units,
the application will not function properly.
