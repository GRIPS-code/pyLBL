from os.path import isfile
from sqlite3 import OperationalError

from numpy import zeros

from pyrad.lbl.continua import OzoneContinuum, WaterVaporForeignContinuum, \
                               WaterVaporSelfContinuum
from pyrad.lbl.hitran import Hitran, Voigt
from pyrad.lbl.hitran.collision_induced_absorption import HitranCIA
from pyrad.lbl.hitran.cross_sections import HitranCrossSection
from pyrad.lbl.tips import TotalPartitionFunction
from pyrad.optics.gas import Gas as Lines


class Gas(object):
    """Gas, capable of calculating its absorption coefficient spectrum.

    Attributes:
        cia: Dictionary where the key is the string chemical formula of the broadener
             and the value is a HitranCIA object.
        continua: List of continuum objects (WaterVaporForeignContinuum and
                  WaterVaporSelfContinuum objects for H2O, OzoneContinuum for O3).
        cross_section: HitranCrossSection object.
        lines: HitranLines object
        molecule: String chemical formula.
    """
    def __init__(self, molecule, broadeners, database):
        """Initializes object, by reading data from the input database.

        Args:
            molecule: String chemical formula.
            broadeners: List of string chemical formulae.
            database: Path to database.
        """
        self.molecule = molecule
        self.lines = Lines(molecule, hitran_database=database, tips_database=database)
        try:
            self.cross_section = HitranCrossSection(molecule, database=database)
        except OperationalError:
            pass
        self.cia = {}
        for broadener in broadeners:
            try:
                self.cia[broadener] = HitranCIA(molecule, broadener, database=database)
            except OperationalError:
                continue
        if molecule == "H2O":
            self.continua = [x(database=database) for x  in [WaterVaporForeignContinuum,
                                                             WaterVaporSelfContinuum]]
        elif molecule == "O3":
            self.continua = [OzoneContinuum(database=database),]
        else:
            self.continua = []

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid):
        """Calculates absorption coefficients.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Dictionary where the keys are string molecular formulae
                                 and the values are volume-mixing ratios [mol mol-1].
            grid: Array of wavenumbers [cm-1].

        Returns:
            Array of absorption coefficients [m-1].
        """
        k = zeros(grid.size)
        k += self.lines.absorption_coefficient(temperature, pressure,
                                               volume_mixing_ratio[self.molecule], grid)
        try:
            k += self.cross_section.absorption_coefficient(temperature, pressure, grid)
        except AttributeError:
            pass
        self.cia = {}
        for broadener, cia in self.cia.items():
            k += cia.absorption_coefficient(temperature, grid) * \
                 number_density(pressure, temperature, volume_mixing_ratio[broadener])
        for continuum in self.continua:
            try:
                k += continuum.absorption_coefficient(temperature, pressure,
                                                      volume_mixing_ratio[self.molecule],
                                                      grid)
            except TypeError:
                k += continuum.absorption_coefficient(grid)
        return k*number_density(pressure, temperature, volume_mixing_ratio[self.molecule])


class Spectroscopy(object):
    """Line-by-line gas optics."""
    def __init__(self, molecules, database=None):
        """Initializes object.  Reads databases to discover what information is available.

        Args:
            local_path: Local directory where database is stored.
        """
        self.database = database
        self.molecules = molecules
        if database is not None and not isfile(self.database):
            self.create_database(molecules, database)

    @staticmethod
    def create_database(molecules, database):
        for molecule in molecules:
            Hitran(molecule, Voigt()).create_database(database)
            TotalPartitionFunction(molecule).create_database(database)
            try:
                HitranCrossSection(molecule).create_database(database)
            except ValueError:
                pass
            for broadener in molecules:
                try:
                    HitranCIA(molecule, broadener).create_database(database)
                except ValueError:
                    continue
            if molecule == "H2O":
                for continuum in [WaterVaporForeignContinuum, WaterVaporSelfContinuum]:
                    continuum().create_database(database)
            elif molecule == "O3":
                OzoneContinuum().create_database(database)

    def load_spectral_inputs(self, wavenumber_lims=None, isotopologues=None):
        """Load spectral data into memory.

        Args:
            wavenumber_lims: Tuple (lower, upper) defining the spectral range [cm-1] for the
                             calculation.  Not used yet.
            isotopologues: Not used yet.
        """
        self.gas = {x: Gas(x, self.molecules, self.database) for x in self.molecules}

    def list_molecules(self):
        """Provides information about each molecule in the local database.

        Returns:
            Dictionary with one entry per gas available in the database
            The value for each gas is a dict with possible (Boolean) keys has_self_continuum, has_foregn_continiuum,
            and collision_induced_absorption (the values are the gas_id(s) of the broadening gas(es))

            I.e.: {"ch4":{},
                   "h2o":{"self_continuum":True,"foregn_continiuum":True},
                   "o2":{"collision_induced_absorption":["o2", "n2"]}
                   }
        """
        return {name: gas.__dict__ for name, gas in self.gas.items()}

    def compute_absorption(self, atmosphere, grid):
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
        k = {x: zeros((atmosphere.pressure.size, grid.size)) for x in self.gas.keys()}
        for name, gas in self.gas.items():
            for i, (temperature, pressure, vmr) in enumerate(atmosphere.points()): 
                k[name][i, :] = gas.absorption_coefficient(temperature, pressure, vmr, grid)
        return k

def number_density(pressure, temperature, volume_mixing_ratio):
    """Calculates the nubmer density using the ideal gas law.

    Args:
        pressure: Pressure [Pa].
        temperature: Temperature [K].
        volume_mixing_ratio: Volume-mixing ratio [mol mol-1].

    Returns:
        Number density [m-3].
    """
    kb = 1.38064852e-23  # Boltzmann constant [J K-1].
    return pressure*volume_mixing_ratio/(kb*temperature)


if __name__ == "__main__":
    from os.path import join

    from netCDF4 import Dataset
    from numpy import arange

    from circ import Atmosphere, Circ

    grid = arange(1., 3250., 1.)
    molecules = ["CH4", "CO", "CO2", "H2O", "N2", "N2O", "O2", "O3"]
    spectroscopy = Spectroscopy(molecules, "test.db")
    spectroscopy.load_spectral_inputs()
    print(spectroscopy.list_molecules())
    circ_data = Circ("circ-case1.nc")
    atmos = Atmosphere(circ_data.temperature, circ_data.pressure, circ_data.vmr)
    absorption_coefficient = spectroscopy.compute_absorption(atmos, grid)
    with Dataset("results.nc", "w") as dataset:
        for name, units, data in zip(["pressure", "wavenumber"], ["Pa", "cm-1"],
                                     [atmos.pressure, grid]):
            dataset.createDimension(name, data.size)
            v = dataset.createVariable(name, "f8", (name,))
            v.setncattr("units", units)
            v[:] = data[:]
        for name, data in absorption_coefficient.items():
            v = dataset.createVariable("{}_absorption_coefficient".format(name), "f8",
                                       ("pressure", "wavenumber"))
            v.setncattr("units", "m-1")
            v[...] = data[...]
