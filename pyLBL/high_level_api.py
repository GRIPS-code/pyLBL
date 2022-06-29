from numpy import sum as npsum
from numpy import unravel_index, zeros
from xarray import DataArray, Dataset

from .atmosphere import Atmosphere
from .database import AliasNotFoundError, CrossSectionNotFoundError, \
                      IsotopologuesNotFoundError, TipsDataNotFoundError, \
                      TransitionsNotFoundError
from .low_level_api import continua, cross_sections, molecular_lines


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


class MoleculeCache(object):
    """Helper class that caches molecular data so it can be reused.

    Attributes:
        gas: Object controlling the current lines backend plugin.
        gas_continua: List of objects controlling the curret continua backend plugin.
    """
    def __init__(self, name, grid, lines_database, lines_engine, continua_engine,
                 cross_sections_engine):
        try:
            self.gas = lines_engine(lines_database, name)
        except (AliasNotFoundError, IsotopologuesNotFoundError,
                TipsDataNotFoundError, TransitionsNotFoundError):
            self.gas = None
        if name == "H2O":
            names = ["{}{}".format(name, x) for x in ["Foreign", "Self"]]
        else:
            names = [name, ]
        try:
            self.gas_continua = [continua_engine[x]() for x in names]
        except KeyError:
            self.gas_continua = None
        try:
            self.cross_section = cross_sections_engine(name, lines_database.arts_crossfit(name))
        except (AliasNotFoundError, CrossSectionNotFoundError):
            self.cross_section = None


class Spectroscopy(object):
    """Line-by-line gas optics.

    Attributes:
        atmosphere: Atmosphere object describing atmospheric conditions.
        cache: Dictionary of MoleculeCache objects.
        continua_backend: String name of model to use for the continua.
        continua_engine: Object exposed by the current molecular continuum backed.
        cross_section_backend: String name of model to use for the cross sections.
        cross_section_engine: Object exposed by the current cross sections backed.
        grid: Numpy array describing the spectral grid [cm-1].
        lines_backend: String name of model to use for lines calculation.
        lines_database: Database object controlling the spectral database.
        lines_engine: Object exposed by the current molecular lines backend.
    """
    def __init__(self, atmosphere, grid, database, mapping=None,
                 lines_backend="pyLBL", continua_backend="mt_ckd",
                 cross_sections_backend="arts_crossfit"):
        """Initializes object.

        Example:
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
            database: Database object to use.
            mapping: Dictionary describing atmospheric dataset variable names.
            lines_backend: String name of model to use for lines calculation.
            continua_backend: String name of model to use for molecular continua.
        """
        self.atmosphere = Atmosphere(atmosphere, mapping=mapping)
        self.grid = grid
        self.lines_database = database
        self.lines_backend = lines_backend
        self.lines_engine = molecular_lines[lines_backend]
        self.continua_backend = continua_backend
        self.continua_engine = continua[continua_backend]
        self.cross_sections_backend = cross_sections_backend
        self.cross_sections_engine = cross_sections[cross_sections_backend]
        self.cache = {}

    def list_molecules(self):
        """Provides a list of molecules available in the specral lines database.

        Returns:
            List of string molecule formulae available in the specral lines database.
        """
        return self.lines_database.molecules()

    def compute_absorption(self, output_format="all", remove_pedestal=None):
        """Computes absorption coefficient [m-1] at specified wavenumbers given temperature,
           pressure, and gas concentrations.

        Args:
            output_format: String describing how the absorption data should be output.
                           "all" - returns absorption spectra from all components (lines,
                                   continuum, cross section) for all gases separately.
                           "gas" - returns total absorption spectra for all gases
                                   separately.
                           "total" - returns the total spectra.
            remove_pedestal: Flag that allows the user to not subtract off the
                             MT-CKD water vapor "pedestal" if desired.

        Returns:
            An xarray Dataset of absorption coefficients [m-1].
        """
        p = self.atmosphere.pressure
        t = self.atmosphere.temperature

        # Initialize the output dataset.
        output = {
            "wavenumber": DataArray(self.grid, dims=("wavenumber",), attrs={"units": "cm-1"}),
        }
        mechanisms = ["lines", "continuum", "cross_section"]
        dims = list(t.dims) + ["mechanism", "wavenumber", ]
        sizes = [x for x in t.sizes.values()] + [len(mechanisms), self.grid.size, ]
        if output_format == "all":
            output["mechanism"] = DataArray(mechanisms, dims=("mechanism",))

        # Calculate the absorption for each molecule at each atmospheric grid point.
        beta = {}
        units = {"units": "m-1"}
        if remove_pedestal is None:
            remove_pedestal = self.continua_backend == "mt_ckd"
        for name, mole_fraction in self.atmosphere.gases.items():
            varname = "{}_absorption".format(name)
            beta[varname] = DataArray(zeros(sizes), dims=dims, attrs=units)
            try:
                # Grab cached spectral database data.
                data = self.cache[name]
            except KeyError:
                # If not already cached, then cache it.
                data = MoleculeCache(name, self.grid, self.lines_database,
                                     self.lines_engine, self.continua_engine,
                                     self.cross_sections_engine)
                self.cache[name] = data
            for i in range(t.data.size):
                vmr = {x: y.data.flat[i] for x, y in self.atmosphere.gases.items()}
                n = number_density(t.data.flat[i], p.data.flat[i],
                                   mole_fraction.data.flat[i])
                j = unravel_index(i, t.data.shape)

                # Calculate lines.
                if data.gas is not None:
                    k = data.gas.absorption_coefficient(t.data.flat[i], p.data.flat[i],
                                                        mole_fraction.data.flat[i], self.grid,
                                                        remove_pedestal=remove_pedestal)
                    indices = tuple(list(j) + [0, slice(None)])
                    beta[varname].values[indices] = n*k[:self.grid.size]

                # Calculate continua.
                if data.gas_continua is not None:
                    indices = tuple(list(j) + [1, slice(None)])
                    for continuum in data.gas_continua:
                        k = continuum.spectra(t.data.flat[i], p.data.flat[i], vmr, self.grid)
                        beta[varname].values[indices] += k[:]

                # Calculate the cross section.
                if data.cross_section is not None:
                    k = data.cross_section.absorption_coefficient(self.grid, t.data.flat[i],
                                                                  p.data.flat[i])
                    indices = tuple(list(j) + [2, slice(None)])
                    beta[varname].values[indices] = n*k[:]

        # Combine the output into the desired form.
        if output_format == "all":
            output.update(beta)
        elif output_format == "gas":
            dims.pop(-2)
            output.update({x: DataArray(npsum(y.values, axis=-2), dims=dims, attrs=units)
                           for x, y in beta.items()})
        else:
            dims.pop(-2)
            data = [DataArray(npsum(x.values, axis=-2), dims=dims, attrs=units) for x in
                    beta.values()]
            output["absorption"] = DataArray(sum([x.values for x in data]),
                                             dims=dims, attrs=units)
        return Dataset(output)
