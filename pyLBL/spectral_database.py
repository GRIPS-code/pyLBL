from logging import info
from os.path import isfile, join

from hapi2 import fetch_isotopologues, fetch_molecules, fetch_parameter_metas, \
                  fetch_transitions, init, Molecule, Transition


_default_config = {
   "engine": "sqlite",
   "database": "local",
   "user": "root",
   "database_dir": ".",
   "host": "http://hitran.org",
   "api_version": "v2",
   "tmpdir": "tmp",
   "api_key": None,
}


def _download_data(molecule, numin, numax, line_list):
    """Downloads isotopologues and transitions for a molecule and stores them in a database.

    Args:
        molecule: HAPI2 Molecule object.
        numin: Wavenumber lower bound (inclusive) [cm-1].
        numax: Wavenumber upper bound (inclusive) [cm-1].
        line_list: Name of temporary files.
    """
    if str(molecule) in ["Chlorine Nitrate",]: return
    try:
        isotopologues = fetch_isotopologues(molecule)
        fetch_transitions(isotopologues, float(numin), float(numax), line_list)
    except Exception as e:
        if str(e) != "Failed to retrieve data for given parameters.": raise


class SpectralDatabase(object):
    """Container the spectral database.

    Attributes:
        database: Path to the local spectral database.
        isotopologues: Dictionary of HAPI2 Isotopologue objects.
        molecules: List of HAPI2 Molecule objects.
        transitions: Dictionary of HAPI2 Transition objects.
    """
    def __init__(self, config, molecules=None, numin=0., numax=60000.,
                 update_database=False):
        """Initialize an object that can interact with the HITRAN database.

        Args:
            config: Dictionary of HAPI configuration settings.
            molecules: List of molecule names.
            numin: Lower bound for transition wavenumbers [cm-1].
            numax: Upper bound for transition wavenumbers [cm-1].
            update_database: Flag that forces the database to be recreated.
        """
        config_ = _default_config.copy()
        config_.update(config)
        self.database = join(config_["database_dir"], "{}.db".format(config_["database"]))
        local_db = isfile(self.database) and not update_database
        init(config_)
        fetch_parameter_metas()
        if local_db:
            info("{} exists, using as local spectral database.".format(self.database))
        else:
            info("Creating local spectral database {}".format(self.database))
            self._create_local_database(molecules, numin, numax)
        self.molecules = Molecule.all()
        self.isotopologues = {}
        self.transitions = {}
        for molecule in self.molecules:
            if str(molecule) in ["Chlorine Nitrate",]: continue
            self.isotopologues[str(molecule)] = Molecule(str(molecule)).isotopologues
            self.transitions[str(molecule)] = Molecule(str(molecule)).transitions

    @staticmethod
    def _create_local_database(molecules, numin, numax, line_list="line-list"):
        """Downloads HITRAN line-by-line data and creates a local SQL database.

        Args:
            molecules: List of molecule names.
            numin: Wavenumber lower bound (inclusive) [cm-1].
            numax: Wavenumber upper bound (inclusive) [cm-1].
            line_list: Name of temporary files.

        Raises:
            ValueError if input molecule is not found in the database.
        """
        all_molecules = fetch_molecules()
        if molecules is None:
            for molecule in all_molecules:
                _download_data(molecule, numin, numax, line_list)
        else:
            for name in molecules:
                for molecule in all_molecules:
                    if name in [str(x) for x in molecule.aliases]:
                        break
                else:
                    raise ValueError("{} not found in spectral database.".format(name))
                _download_data(molecule, numin, numax, line_list)


    @staticmethod
    def load_line_parameters(formula, numin, numax):
        """Reads the HITRAN molecular line parameters from a local SQL database.

        Args:
            formula: String chemical formula (i.e. H2O).
            numin: Wavenumber lower bound (inclusive) [cm-1].
            numax: Wavenumber upper bound (inclusive) [cm-1].

        Returns:
            A list of Transition objects.
        """
        return Molecule(formula).transitions. \
               filter(Transition.nu>=numin).filter(Transition.nu<=numax)