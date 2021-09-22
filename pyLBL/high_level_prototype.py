from xarray import DataArray, Dataset

from hapi2 import Molecule
from numpy import unravel_index, zeros

from .atmosphere import Atmosphere
from .low_level_prototype import continua, molecular_lines
from .spectral_database import Hapi


kb = 1.38064852e-23  # Boltzmann constant [J K-1].


def number_density(temperature, pressure, volume_mixing_ratio):
    """Calculates the number density using the ideal gas law.

    Args:
        temperature: Temperature [K].
        pressure: Pressure [Pa].
        volume_mixing_ratio: Volume-mixing ratio [mol mol-1].

    Returns:
        Number density [m-3].
    """
    return pressure*volume_mixing_ratio/(kb*temperature)


class Spectroscopy(object):
    """Line-by-line gas optics."""
    def __init__(self, lines_backend="grtcode", continua_backend="mt_ckd"):
        """Initializes object.  Reads databases to discover what information is available.

        Args:
            local_path: Local directory where database is stored.
        """
        self.lines_database = Hapi()
        self.lines_engine = molecular_lines[lines_backend]
        self.continua_engine = continua[continua_backend]

#   def load_spectral_inputs(self, atmosphere, spectral_grid):
#       """Load spectral data into memory.

#       Args:
#           wavenumber_lims: Tuple (lower, upper) defining the spectral range [cm-1] for the
#                            calculation.  Not used yet.
#           isotopologues: Not used yet.
#       """
#       for gas in Atmosphere(atmosphere).keys():
#           self.lines_database.load_line_parameters(gas, spectral_grid[0], spectral_grid[-1])

    def list_molecules(self):
        """Provides information about each molecule in the local database.

        Returns:
            Dictionary where the key is name of each gas available in the database
            and the value is a dictionary describing the different objects capable
            of calculating absorption coefficients.
        """
        return [str(x) for x in self.lines_database.molecules]

    def compute_absorption(self, atmosphere, grid, mapping=None):
        """Computes absorption coefficient [m-1] at specified wavenumbers given temperature,
           pressure, and gas concentrations 

        Args:
            atmosphere: an xarray Dataset.
            grid: Wavenumber grid array [cm-1].
            mapping: Dictionary describing dataset variable names.  Should have the form:
                     {"play": <name of pressure variable in dataset>,
                      "tlay": <name of temperature variable in dataset>,
                      "mole_fraction:
                          {"H2O" : <name of water vapor mole fraction variable in dataset>,
                           "CO2" : <name of carbon dioxided mole fraction variable in dataset>,
                           ...
                          }
                     }

        Returns:
            An xarray Dataset of absorption coefficients [m-1].
        """
        atm = Atmosphere(atmosphere, mapping=mapping)
        p = atm.pressure
        t = atm.temperature
        dims = list(t.dims) + ["mechanism", "wavenumber"]
        sizes = tuple([x for x in t.sizes.values()] + [2, grid.size])
        beta = {"wavenumber": DataArray(grid, dims=("wavenumber",), attrs={"units": "cm-1"}),
                "mechanism": DataArray(["lines", "continuum"], dims=("mechanism",))}
        for name, mole_fraction in atm.gases.items():
            varname = "{}_absorption".format(name)
            beta[varname] = DataArray(zeros(sizes), dims=dims, attrs={"units": "m-1"})
            avg_mass = sum([x.abundance*x.mass for x in Molecule(name).isotopologues])
            mol_id = Molecule(name).id
            num_iso = len(Molecule(name).isotopologues)
            transitions = self.lines_database.load_line_parameters(name, grid[0], grid[-1])
            gas = self.lines_engine(transitions, mol_id, num_iso, avg_mass)

            names = ["{}{}".format(name, x) for x in ["Foreign", "Self"]] \
                    if name == "H2O" else [name,]
            try:
                gas_continua = [self.continua_engine[x]() for x in names]
            except KeyError:
                gas_continua = None
            for i in range(t.data.size):
                vmr = {x: y.data.flat[i] for x, y in atm.gases.items()}
                n = number_density(t.data.flat[i], p.data.flat[i],
                                   mole_fraction.data.flat[i])
                j = unravel_index(i, t.data.shape)

                # Calculate lines.
                k = gas.absorption_coefficient(t.data.flat[i], p.data.flat[i],
                                               mole_fraction.data.flat[i], grid)
                indices = tuple(list(j) + [0, slice(None)])
                beta[varname].values[indices] = n*k[:]

                # Calculate continua.
                if gas_continua is not None:
                    indices = tuple(list(j) + [1, slice(None)])
                    for continuum in gas_continua:
                        k = continuum.spectra(t.data.flat[i], p.data.flat[i], vmr, grid)
                        beta[varname].values[indices] += k[:]
        return Dataset(beta)
