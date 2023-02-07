from logging import getLogger
from os.path import join
from pathlib import Path
from subprocess import run
from sys import stderr, stdout

from numpy import asarray


info = getLogger(__name__).info
arts_installed = False
try:
    from pyarts.workspace import Workspace
    arts_installed = True
except ImportError:
    info("pyarts is not installed.")


def download_data(cwd=None):
    """Downloads molecular line data.

    Args:
        cwd: Directory to run the download in.
    """
    base_url = "https://arts.mi.uni-hamburg.de/svn/rt/arts-xml-data/trunk"
    for url in ["spectroscopy/Hitran/", ]:
        info("Downloading HITRAN data from {}/{}.".format(base_url, url))
        run(["svn", "co", "-q", "/".join([base_url, url])], stdout=stdout,
            stderr=stderr, check=True, cwd=cwd)


def load_data():
    """Downloads molecular line data if not found in the package directory.

    Returns:
        Absolute path of the directory containing molecular line data.
    """
    pkg_dir = Path(__file__).parent
    hitran = pkg_dir / "Hitran"
    if not (hitran.exists() and hitran.is_dir()):
        download_data(cwd=str(pkg_dir))
    return str(hitran.absolute())


def configure_workspace(verbosity=0):
    """Configures the ARTS application.

    Args:
        verbosity: ARTS verbosity level.

    Returns:
        A Workspace object.
    """
    workspace = Workspace(verbosity=0)
    for name in ["general", "continua", "agendas"]:
        workspace.execute_controlfile(join("general", "{}.arts".format(name)))
    workspace.verbositySetScreen(workspace.verbosity, verbosity)
    workspace.jacobianOff()
    workspace.Copy(workspace.abs_xsec_agenda, workspace.abs_xsec_agenda__noCIA)
    workspace.AtmosphereSet1D()
    return workspace


class PyArtsGas(object):
    def __init__(self, lines_database, formula):
        if not arts_installed:
            raise ValueError("pyarts is not installed.")
        hitran_directory = "{}/".format(load_data())
        self.formula = formula
        self.workspace = configure_workspace(verbosity=2)
        self.workspace.abs_speciesSet(species=[formula])
        self.workspace.ArrayOfIndexSet(self.workspace.abs_species_active, [0])
        self.workspace.abs_lines_per_speciesReadSpeciesSplitCatalog(basename=hitran_directory)
        self.workspace.abs_lines_per_speciesSetCutoff(option="ByLine", value=750.e9)
        self.workspace.ArrayOfArrayOfAbsorptionLinesCreate("abs_lines_per_species_backup")
        self.workspace.Copy(self.workspace.abs_lines_per_species_backup,
                            self.workspace.abs_lines_per_species)
        self.workspace.isotopologue_ratiosInitFromBuiltin()

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid,
                               remove_pedestal=False, cut_off=25):
        """Calculates absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Numpy array defining the spectral grid [cm-1].
            remove_pedestal: Flag specifying if a pedestal should be subtracted.
            cut_off: Wavenumber cut-off distance [cm-1] from line centers.

        Returns:
            Numpy array of absorption coefficients [m2].
        """
        # Configure spectral grid.
        self.workspace.f_grid = grid
        self.workspace.FrequencyFromCGSKayserWavenumber(self.workspace.f_grid,
                                                        self.workspace.f_grid)
        self.workspace.abs_lines_per_speciesCompact()

        # Configure the atmosphere.
        self.workspace.NumericSet(self.workspace.rtp_pressure, pressure)
        self.workspace.NumericSet(self.workspace.rtp_temperature, temperature)
        self.workspace.VectorSet(self.workspace.rtp_vmr, asarray([volume_mixing_ratio]))
        self.workspace.Touch(self.workspace.abs_nlte)
        self.workspace.AbsInputFromRteScalars()

        # Calculate the absorption coefficient.
        self.workspace.lbl_checkedCalc()
        self.workspace.abs_xsec_agenda_checkedCalc()
        self.workspace.abs_xsec_per_speciesInit()
        self.workspace.abs_xsec_per_speciesAddLines()
        self.workspace.Copy(self.workspace.abs_lines_per_species,
                            self.workspace.abs_lines_per_species_backup)
        return self.workspace.abs_xsec_per_species.value.copy()[0][:, 0]
