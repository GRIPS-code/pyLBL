Installation
------------

Requirements
~~~~~~~~~~~~

:code:`pyLBL` requires python 3 (>= version 3.6) and requires the following:

* matplotlib
* netCDF4
* numpy
* scipy
* sqlalchemy
* xarray

Installing from Conda Forge
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The preferred way to install :code:`pyLBL` is with conda_.  Users may it find it useful
to create their own conda environment and install the model inside it by running:

.. code-block:: bash

  # Optional: If the user would like to set up a custom conda environment.
  conda config --add envs_dirs <location>
  conda config --add pkgs_dirs <location>
  conda create -n <name>
  conda activate <name>

  conda install -c conda-forge pyLBL

Here the :code:`<name>` of the environment can be anything the user desires and the
:code:`<location>` paths can be in user-writeable locations (if for example users are
not allowed to write to system directories).

Installing from Github
~~~~~~~~~~~~~~~~~~~~~~

:code:`pyLBL` can also be obtained from github by running:

.. code-block:: bash

  git clone --recursive https://github.com/GRIPS-code/pyLBL.git
  cd pyLBL
  python3 setup.py install

In order to contribute, please fork the repository and submit issues and pull requests.

.. _conda: https://anaconda.org/conda-forge/pylbl
