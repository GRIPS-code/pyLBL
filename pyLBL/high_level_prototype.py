from xarray import DataArray, Dataset

from hapi2 import Molecule
from numpy import unravel_index, zeros

from .atmosphere import Atmosphere
from .low_level_prototype import continua, molecular_lines
from .spectral_database import SpectralDatabase


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
    """Line-by-line gas optics.

    Attributes:
        atmosphere: Atmosphere object describing atmospheric conditions.
        grid: Numpy array describing the spectral grid [cm-1].
        lines_database: SpectralDatabase object controlling the spectral database.
        lines_engine: Object exposed by the current molecular lines backend.
        continua_engine: Object exposed by the current molecular cotinuum backed.
    """
    def __init__(self, atmosphere, grid, database=None, mapping=None, hapi_config=None,
                 lines_backend="grtcode", continua_backend="mt_ckd"):
        """Initializes object.

        Example::
            mapping = {
                "play": <name of pressure variable in dataset>,
                "tlay": <name of temperature variable in dataset>,
                "mole_fraction: {
                    "H2O" : <name of water vapor mole fraction variable in dataset>,
                    "CO2" : <name of carbon dioxided mole fraction variable in dataset>,
                    ...
                }
            }

        Args:
            atmosphere: xarray Dataset describing atmospheric conditions.
            grid: Wavenumber grid array [cm-1].
            database: SpectralDatabase object to use.
            mapping: Dictionary describing atmospheric dataset variable names.
            hapi_config: Dictionary of HAPI2 configuration options.
            lines_backend: String name of model to use for lines calculation.
            continua_backend: String name of model to use for molecular continua.
        """
        self.atmosphere = Atmosphere(atmosphere, mapping=mapping)
        self.grid = grid
        self.lines_database = database if database is not None else \
            SpectralDatabase(molecules=self.atmosphere.gases.keys(),
                             numin=grid[0], numax=grid[-1], config=hapi_config)
        self.lines_engine = molecular_lines[lines_backend]
        self.continua_engine = continua[continua_backend]

    def list_molecules(self):
        """Provides information about each molecule in the local database.

        Returns:
            Dictionary where the key is name of each gas available in the database
            and the value is a dictionary describing the different objects capable
            of calculating absorption coefficients.
        """
        return [str(x) for x in self.lines_database.molecules]

    def compute_absorption(self):
        """Computes absorption coefficient [m-1] at specified wavenumbers given temperature,
           pressure, and gas concentrations 

        Returns:
            An xarray Dataset of absorption coefficients [m-1].
        """
        p = self.atmosphere.pressure
        t = self.atmosphere.temperature
        dims = list(t.dims) + ["mechanism", "wavenumber"]
        sizes = tuple([x for x in t.sizes.values()] + [2, self.grid.size])
        beta = {"wavenumber": DataArray(self.grid, dims=("wavenumber",), attrs={"units": "cm-1"}),
                "mechanism": DataArray(["lines", "continuum"], dims=("mechanism",))}
        for name, mole_fraction in self.atmosphere.gases.items():
            varname = "{}_absorption".format(name)
            beta[varname] = DataArray(zeros(sizes), dims=dims, attrs={"units": "m-1"})
            avg_mass = sum([x.abundance*x.mass for x in Molecule(name).isotopologues])
            mol_id = Molecule(name).id
            num_iso = len(Molecule(name).isotopologues)
            transitions = self.lines_database.load_line_parameters(name, self.grid[0], self.grid[-1])
            gas = self.lines_engine(transitions, mol_id, num_iso, avg_mass)

            names = ["{}{}".format(name, x) for x in ["Foreign", "Self"]] \
                    if name == "H2O" else [name,]
            try:
                gas_continua = [self.continua_engine[x]() for x in names]
            except KeyError:
                gas_continua = None
            for i in range(t.data.size):
                vmr = {x: y.data.flat[i] for x, y in self.atmosphere.gases.items()}
                n = number_density(t.data.flat[i], p.data.flat[i],
                                   mole_fraction.data.flat[i])
                j = unravel_index(i, t.data.shape)

                # Calculate lines.
                k = gas.absorption_coefficient(t.data.flat[i], p.data.flat[i],
                                               mole_fraction.data.flat[i], self.grid)
                indices = tuple(list(j) + [0, slice(None)])
                beta[varname].values[indices] = n*k[:]

                # Calculate continua.
                if gas_continua is not None:
                    indices = tuple(list(j) + [1, slice(None)])
                    for continuum in gas_continua:
                        k = continuum.spectra(t.data.flat[i], p.data.flat[i], vmr, self.grid)
                        beta[varname].values[indices] += k[:]
        return Dataset(beta)
