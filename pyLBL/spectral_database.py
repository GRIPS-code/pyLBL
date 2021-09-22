from logging import info
from os.path import isfile

from hapi2 import fetch_isotopologues, fetch_molecules, fetch_parameter_metas, \
                  fetch_transitions, init, Molecule, Transition


default_config = {
   "engine": "sqlite",
   "database": "local",
   "user": "root",
   "database_dir": ".",
   "host": "http://hitran.org",
   "api_version": "v2",
   "tmpdir": "tmp",
   "api_key": None,
}


class Hapi(object):
    def __init__(self, config, update_database=False):
        """Initialize an object that can interact with the HITRAN database.

        Args:
            config: Dictionary of HAPI configuration settings.
            update_database: Flag that forces the database to be recreated.
        """
        hapi_config = default_config.copy()
        hapi_config.update(config)
        init(hapi_config)
        self.database = "{}.db".format(hapi_config["database"])
        fetch_parameter_metas()
        if isfile(self.database) and not update_database:
            info("{} exists, using as local spectral database.".format(self.database))
        else:
            info("Creating local spectral database {}".format(self.database))
            self._create_local_database()
        self.molecules = Molecule.all()
        self.isotopologues = {}
        self.transitions = {}
        for molecule in self.molecules:
            if str(molecule) in ["Chlorine Nitrate",]: continue
            self.isotopologues[str(molecule)] = Molecule(str(molecule)).isotopologues
            self.transitions[str(molecule)] = Molecule(str(molecule)).transitions

    @staticmethod
    def _create_local_database(numin=0, numax=60000, line_list="line-list"):
        """Downloads HITRAN line-by-line data and creates a local SQL database.

        Args:
            numin: Wavenumber lower bound (inclusive) [cm-1].
            numax: Wavenumber upper bound (inclusive) [cm-1].
            line_list: Name of temporary files.
        """
        molecules = fetch_molecules()
        for molecule in molecules:
            if str(molecule) in ["Chlorine Nitrate",]: continue
            try:
                isotopologues = fetch_isotopologues(molecule)
                fetch_transitions(isotopologues, numin, numax, line_list)
            except Exception as e:
                if str(e) != "Failed to retrieve data for given parameters.": raise

    @staticmethod
    def load_line_parameters(formula, numin, numax):
        """Reads the HITRAN molecular line parameters from a local SQL database.

        Args:
            formula: String chemical formula (i.e. H2O).

        Returns:
            A list of Transition objects.
        """
        return Molecule(formula).transitions. \
               filter(Transition.nu>=numin).filter(Transition.nu<=numax)
