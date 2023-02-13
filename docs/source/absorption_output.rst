Absorption Output
-----------------

Absorption coefficients can be calculated using the :code:`Spectroscopy` object (described
in the previous API section) by running:

.. code-block:: python

  # Calculate the absorption coefficients [m-1] and return them in an xarray Dataset.
  absorption = spectroscopy.compute_absorption(output_format="all")

  # Optional: convert dataset to netcdf.
  absorption.to_netcdf("<name of output file>")

The output is returned as an xarray :code:`Dataset`.  The absorption coefficients are
calculated in units of inverse meters, so that optical depths can be calculated by the user
by integrationg them over any desired path.  The exact format of the output data
depends on the value of the :code:`output_format` argument.  When set to :code:`"all"` (which
is currently the default), the dataset will return the spectra split up by molecule
and mechansim (lines, continuum, cross_section). If a mechanism does not apply to a specific
molecule, then values are set to zero accordingly.  An example viewed in netCDF format
(i.e., like by running :code:`ncdump` on the dataset) would look like this:

.. code-block:: none

  netcdf absorption {
  dimensions:
          wavenumber = 49990 ;
          mechanism = 3 ;
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
          mechanism = "lines", "continuum", "cross_section" ;
  }

If the :code:`output_format` argument is instead set to :code:`"gas"`, the spectra for
the different mechanims will be summed for each molecule, yielding output that looks
like this (in netCDF format):

.. code-block:: none

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

Lastly, if the :code:`output_format` argument is set to any other value, only the total
absorption spectra (summed over all molecules) will be returned.  In netCDF format, the
resulting dataset will appear like this:

.. code-block:: none

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

We hope that formatted in this way will ease interaction between this package and other
pangeo tools.
