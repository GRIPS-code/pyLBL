from xarray import DataArray, Dataset

from numpy import sum as npsum
from numpy import unravel_index, zeros

from .atmosphere import Atmosphere
from .low_level_api import continua, molecular_lines


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
    def __init__(self, name, grid, lines_database, lines_engine, continua_engine):
        self.gas = lines_engine(*lines_database.gas(name))
        names = ["{}{}".format(name, x) for x in ["Foreign", "Self"]] \
                 if name == "H2O" else [name,]
        try:
            self.gas_continua = [continua_engine[x]() for x in names]
        except KeyError:
            self.gas_continua = None


class Spectroscopy(object):
    """Line-by-line gas optics.

    Attributes:
        atmosphere: Atmosphere object describing atmospheric conditions.
        cache: Dictionary of MoleculeCache objects.
        continua_engine: Object exposed by the current molecular cotinuum backed.
        grid: Numpy array describing the spectral grid [cm-1].
        lines_database: Database object controlling the spectral database.
        lines_engine: Object exposed by the current molecular lines backend.
    """
    def __init__(self, atmosphere, grid, database, mapping=None, 
                 lines_backend="pyLBL", continua_backend="mt_ckd"):
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
            database: Database object to use.
            mapping: Dictionary describing atmospheric dataset variable names.
            lines_backend: String name of model to use for lines calculation.
            continua_backend: String name of model to use for molecular continua.
        """
        self.atmosphere = Atmosphere(atmosphere, mapping=mapping)
        self.grid = grid
        self.lines_database = database
        self.lines_engine = molecular_lines[lines_backend]
        self.continua_engine = continua[continua_backend]
        self.cache = {}

    def list_molecules(self):
        """Provides a list of molecules available in the specral lines database.

        Returns:
            List of string molecule formulae available in the specral lines database.
        """
        return self.lines_database.molecules()

    def compute_absorption(self, output_format="all"):
        """Computes absorption coefficient [m-1] at specified wavenumbers given temperature,
           pressure, and gas concentrations 

        Args:
            output_format: String describing how the absorption data should be output.
                           "all" - returns absorption spectra from all components (lines,
                                   continuum, cross section) for all gases separately.
                           "gas" - returns total absorption spectra for all gases
                                   separately.
                           "total" - returns the total spectra.

        Returns:
            An xarray Dataset of absorption coefficients [m-1].
        """
        p = self.atmosphere.pressure
        t = self.atmosphere.temperature

        # Initialize the output dataset.
        output = {
            "wavenumber": DataArray(self.grid, dims=("wavenumber",), attrs={"units": "cm-1"}),
        }
        mechanisms = ["lines", "continuum"]
        dims = list(t.dims) + ["mechanism", "wavenumber",]
        sizes = [x for x in t.sizes.values()] + [len(mechanisms), self.grid.size,]
        if output_format == "all":
            output["mechanism"] = DataArray(["lines", "continuum"], dims=("mechanism",))

        # Calculate the absorption for each molecule at each atmospheric grid point.
        beta = {}
        units = {"units": "m-1"}
        for name, mole_fraction in self.atmosphere.gases.items():
            varname = "{}_absorption".format(name)
            beta[varname] = DataArray(zeros(sizes), dims=dims, attrs=units)
            try:
                # Grab cached spectral database data.
                data = self.cache[name]
            except KeyError:
                # If not already cached, then cache it.
                data = MoleculeCache(name, self.grid, self.lines_database,
                                     self.lines_engine, self.continua_engine)
                self.cache[name] = data
            for i in range(t.data.size):
                vmr = {x: y.data.flat[i] for x, y in self.atmosphere.gases.items()}
                n = number_density(t.data.flat[i], p.data.flat[i],
                                   mole_fraction.data.flat[i])
                j = unravel_index(i, t.data.shape)

                # Calculate lines.
                k = data.gas.absorption_coefficient(t.data.flat[i], p.data.flat[i],
                                                    mole_fraction.data.flat[i], self.grid)
                indices = tuple(list(j) + [0, slice(None)])
                beta[varname].values[indices] = n*k[:]

                # Calculate continua.
                if data.gas_continua is not None:
                    indices = tuple(list(j) + [1, slice(None)])
                    for continuum in data.gas_continua:
                        k = continuum.spectra(t.data.flat[i], p.data.flat[i], vmr, self.grid)
                        beta[varname].values[indices] += k[:]

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
