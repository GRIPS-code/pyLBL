from os.path import isfile
from sqlite3 import OperationalError

from numpy import log, zeros

from pyrad.lbl.continua import OzoneContinuum, WaterVaporForeignContinuum, \
                               WaterVaporSelfContinuum
from pyrad.lbl.hitran import Hitran, Voigt
from pyrad.lbl.hitran.collision_induced_absorption import HitranCIA
from pyrad.lbl.hitran.cross_sections import HitranCrossSection
from pyrad.lbl.tips import TotalPartitionFunction
from pyrad.optics.gas import Gas as Lines


air_molar_mass = 0.0289647  # Molar mass of air [kg mol-1].
Na = 6.02214076e23  # Avogadro's number [mol-1].
kb = 1.38064852e-23  # Boltzmann constant [J K-1].
g = 9.80  # Acceleration due to gravity [m s-2].


class Absorption(object):
    """Container for absorption coefficients.

    Attributes:
        cia: Dictionary of arrays of absorption coefficients [m-1] from collision-induced
             absorption tables.
        continua: Array of absorption coefficients [m-1] from continua.
        cross_section: Array of absorption coefficients [m-1] from cross section tables.
        lines: Array of absorption coefficients [m-1] from spectral lines.
    """
    def __init__(self, lines, cross_section, cia, continua):
        """Initializes object.

        Args:
            lines: Array of absorption coefficients [m-1] from spectral lines.
            cross_section: Array of absorption coefficients [m-1] from cross section tables.
            cia: Dictionary of arrays of absorption coefficients [m-1] from collision-induced
                 absorption tables.
            continua: Array of absorption coefficients [m-1] from continua.
        """
        self.lines = lines
        self.cross_section = cross_section
        self.cia = cia
        self.continua = continua

    @property
    def total(self):
        """Adds the absorption coefficient components together.

        Returns:
            Array of total absorption coefficients [m-1].
        """
        cia = zeros(self.lines.size)
        for data in self.cia.values():
            cia += data
        return self.lines + self.cross_section + self.continua + cia


class Atmosphere(object):
    """Atmosphere object.

    Attributes:
        temperature: Array of temperatures [K].
        pressure: Array of pressures [Pa].
        volume_mixing_ratio: Dictionary where the keys are the molecule formulae and the
                             values are arrays of volume-mixing ratios [mol mol-1].
    """
    def __init__(self, temperature, pressure, volume_mixing_ratio, z=None):
        """Initializes object.

        Args:
            temperature: Array of temperatures [K].
            pressure: Array of pressures [Pa].
            volume_mixing_ratio: Dictionary where the keys are the molecule formulae and the
                                 values are arrays of volume-mixing ratios [mol mol-1].
        """
        self.temperature = temperature
        self.pressure = pressure
        if z is None:
            # Calculate distance between points assuming the atmosphere is 1D using the
            # hydrostatic approximaxtion and ideal gas law.
            self.z = zeros(temperature.size)
            for i in range(1, self.z.size):
               self.z[i] = self.z[i-1] + hydrostatic_distance(pressure[i-1:i+1],
                                                              temperature[i])
        else:
            self.z = z
        self.volume_mixing_ratio = volume_mixing_ratio

    def points(self):
        for i in range(self.temperature.size):
            yield self.temperature[i], self.pressure[i], \
                  {x: y[i] for x, y in self.volume_mixing_ratio.items()}


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
            Absorption object.
        """
        vmr = volume_mixing_ratio[self.molecule]
        n = number_density(pressure, temperature, vmr)
        k_lines = self.lines.absorption_coefficient(temperature, pressure, vmr, grid)*n
        try:
            k_cross_section = self.cross_section.absorption_coefficient(temperature, pressure,
                                                                        grid)*n
        except AttributeError:
            k_cross_section = zeros(grid.size)
        k_cia = {x: zeros(grid.size) for x in self.cia.keys()}
        for broadener, cia in self.cia.items():
            k_cia[broadener] += cia.absorption_coefficient(temperature, grid)*n * \
                                number_density(pressure, temperature,
                                volume_mixing_ratio[broadener])
        k_continua = zeros(grid.size)
        for continuum in self.continua:
            try:
                k_continua += continuum.absorption_coefficient(temperature, pressure,
                                                               vmr, grid)*n
            except TypeError:
                k_continua += continuum.absorption_coefficient(grid)*n
        return Absorption(k_lines, k_cross_section, k_cia, k_continua)


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
            Dictionary where the key is name of each gas available in the database
            and the value is a dictionary describing the different objects capable
            of calculating absorption coefficients.
        """
        return {name: gas.__dict__ for name, gas in self.gas.items()}

    def compute_absorption(self, atmosphere, grid):
        """Computes absorption coefficient [m-1] at specified wavenumbers given temperture,
           pressure, and gas concentrations 

        Args:
            atmosphere: an Atmosphere object.
            grid: Wavenumber grid array [cm-1].

        Returns:
            An array of absorption coefficients.
        """
        k = {x: zeros((atmosphere.pressure.size, grid.size)) for x in self.gas.keys()}
        for name, gas in self.gas.items():
            for i, (temperature, pressure, vmr) in enumerate(atmosphere.points()): 
                k[name][i, :] = gas.absorption_coefficient(temperature, pressure, vmr,
                                                           grid).total
        return k


def hydrostatic_distance(pressure, temperature, molar_mass=air_molar_mass):
    """Calculates the distance between two pressures assuming hydrostatic equilibrium
       and constant temperature, using the ideal gas law.

    Args:
        pressure: 2-element array of pressures [Pa].
        temperature: Temperature [K].
        molar_mass: Molar mass [kg mol-1].

    Returns:
        Distance [m] between the two atmospheric pressures.
    """
    return (Na*kb*temperature/(molar_mass*g))*log(pressure[1]/pressure[0])


def number_density(pressure, temperature, volume_mixing_ratio):
    """Calculates the number density using the ideal gas law.

    Args:
        pressure: Pressure [Pa].
        temperature: Temperature [K].
        volume_mixing_ratio: Volume-mixing ratio [mol mol-1].

    Returns:
        Number density [m-3].
    """
    return pressure*volume_mixing_ratio/(kb*temperature)
