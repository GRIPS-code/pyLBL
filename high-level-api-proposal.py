class Spectroscopy(object)
    """Line-by-line gas optics."""
    def __init__(self, local_path=None, **kwargs):
        """Initializes object - reads databases to discover what information is available

        Args:
            local_path: Local directory where databases are stored
        """

    def create(self, molecules, wavenumber_lims=None, isotopologues, local_path=None, **kwargs):
        """Downloads or otherwise caches data for the specified molecules

        Args:
            molecules: List of strings describing chemical formula (i.e. "H2O").
            wavenumber_lims: tuple of lowest/highest wavenomber for which data is to be downloaded
            isotopologues: Lists of Isotopologue objects.
            local_path: Local directory where databases are to be stored
        """

        # Call __init__ after downloading data?

    def list_molecues(self, local_path=None, **kwargs):
        """Provides information about each molecule in the local database

        Args:
            local_path: Local directory where databases are stored

        Returns:
            Dictionary with one entry per gas available in the database
            The value for each gas is a dict with possible (Boolean) keys has_self_continuum, has_foregn_continiuum,
            and collision_induced_absorption (the values are the gas_id(s) of the broadening gas(es))

            I.e.: {"ch4":{},
                   "h2o":{"self_continuum":True,"foregn_continiuum":True},
                   "o2":{"collision_induced_absorption":["o2", "n2"]}
                   }
        """

    def compute_absorption(self, molecules, atmosphere, grid,
                           local_path=None, provide_components=False, **kwargs):
        """Computes absorption coefficient (inverse meters per molecule) at specified
           wavenumbers given temperture, pressure, and gas concentrations 

        Args:
            molecules: one or more strings describing chemical formula (i.e. "H2O").
            atmosphere: an object describing temperatures, pressures, and gas concentrations
              One option is an xarray Dataset where the names and/or attributes of the DataArrays specify which variables
              to use, units, etc.
            grid: Wavenumber grid array [cm-1].
            isotopologues: Lists of Isotopologue objects.
            local_path: Local directory where databases are to be stored


        Returns:
            xarray Datasets with absorption coefficient (inverse meters per molecule) for each values of molecule
            on the spectral grid and all other coordinates in the atmosphere objects
            If provide_components we break out the components of absorption: more coordinates? Dicts?
        """
