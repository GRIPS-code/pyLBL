API
---

This application aims to improve on existing line-by-line radiative transfer models
by separating the data management and calculation.  Data management is handled by
a :code:`Database` object, as described in the previous section, to eliminate the need to
explicitly interact with the ascii/data files typically needed when using existing
line-by-line models.  Absorption spectra calculation is handled by
a :code:`Spectroscopy` object, which allows the user to specify which molecular
lines, continua, and cross section models they would like to use.

Absorption calculation
~~~~~~~~~~~~~~~~~~~~~~

A :code:`Spectroscopy` object allow users to choose which models are used to calculate the
molecular lines, various molecular continua, and absorption cross sections.  Currently,
the supported models are as follows:

============== ==========================
Component      Models
============== ==========================
lines          pyLBL, pyarts(in progress)
continua       MT-CKD
cross sections arts-crossfit
============== ==========================

For example, to create a :code:`Spectroscopy` object using the default spectral
lines model and the MT-CKD continuum, use:

.. code-block:: python

  from pyLBL import Spectroscopy

  spectroscopy = Spectroscopy(atmosphere, grid, database, mapping=mapping,
                              lines_backend="pyLBL", continua_backend="mt_ckd",
                              cross_sections_backend="arts_crossfit")

Here the :code:`atmosphere` and :code:`mapping` arguments are xarray :code:`Dataset` and
python dictionary objects respectively that describe the input atmospheric
conditions (see the previous "Atmospheric Inputs" section).
The :code:`database` and :code:`grid` arguments are :code:`Database` and numpy array
objects that describe the spectral inputs (see the previous "Spectral inputs" section).
The last three arguments (:code:`lines_backend`, :code:`continua_backend`,
and :code:`cross_sections_backend`) are strings that determine which models will be
used for each component of the calculate (see the table above).

Absorption coefficients can be calculated using the :code:`Spectroscopy` object
by running:

.. code-block:: python

  absorption = spectroscopy.compute_absorption(output_format="all")

See the next section "Absorption Output" which discusses the output options.
