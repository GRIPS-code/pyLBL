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

  conda create -n env
  conda activate env
  conda install -c conda-forge pyLBL

Here the name of the environment :code:`env` can be anything the user desires.

Installing from Github
~~~~~~~~~~~~~~~~~~~~~~

:code:`pyLBL` can also be obtained from github by running:

.. code-block:: bash

  git clone --recursive https://github.com/GRIPS-code/pyLBL.git
  cd pyLBL
  python3 setup.py install

In order to contribute, please fork the repository and submit issues and pull requests.

.. _conda: https://anaconda.org/conda-forge/pylbl
