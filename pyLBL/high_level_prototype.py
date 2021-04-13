from os.path import isfile
from sqlite3 import OperationalError

from numpy import log, unravel_index, zeros
from xarray import DataArray, Dataset

from pyrad.lbl.continua import OzoneContinuum, WaterVaporForeignContinuum, \
                               WaterVaporSelfContinuum
from pyrad.lbl.hitran import Hitran, Voigt
from pyrad.lbl.hitran.collision_induced_absorption import HitranCIA
from pyrad.lbl.hitran.cross_sections import HitranCrossSection
from pyrad.lbl.tips import TotalPartitionFunction
from pyrad.optics.gas import Gas as Lines


kb = 1.38064852e-23  # Boltzmann constant [J K-1].


class Gas(object):
    """Gas, capable of calculating its absorption coefficient spectrum.

    Attributes:
        cia: Dictionary where the key is the string chemical formula of the broadener
             and the value is a HitranCIA object.
        continua: Dictionary where the key is the continuum type and the value is a
                  continuum object (WaterVaporForeignContinuum and WaterVaporSelfContinuum
                  objects for H2O, OzoneContinuum for O3).
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
        try:
            self.lines = Lines(molecule, hitran_database=database, tips_database=database)
        except OperationalError:
            pass
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
        self.continua = {}
        if molecule == "H2O":
            self.continua["foregin"] = WaterVaporForeignContinuum(database=database)
            self.continua["self"] = WaterVaporSelfContinuum(database=database)
        elif molecule == "O3":
            self.continua["all"] = OzoneContinuum(database=database)

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid):
        """Calculates absorption coefficients.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Dictionary where the keys are string molecular formulae
                                 and the values are volume-mixing ratios [mol mol-1].
            grid: Array of wavenumbers [cm-1].

        Returns:
            Dictionary of absorption coefficients.
        """
        vmr = volume_mixing_ratio[self.molecule]
        n = number_density(pressure, temperature, vmr)
        k = {}
        k["lines"] = self.lines.absorption_coefficient(temperature, pressure, vmr, grid)*n
        try:
            k["cross section"] = self.cross_section.absorption_coefficient(temperature, pressure,
                                                                           grid)*n
        except AttributeError:
            k["cross section"] = zeros(grid.size)
        k["cia"] = {x: zeros(grid.size) for x in self.cia.keys()}
        for broadener, cia in self.cia.items():
            k["cia"][broadener] += cia.absorption_coefficient(temperature, grid)*n * \
                                   number_density(pressure, temperature,
                                                  volume_mixing_ratio[broadener])
        k["continua"] = {x: zeros(grid.size) for x in self.continua.keys()}
        for mechanism, continuum in self.continua.items():
            try:
                k["continua"][mechanism] += continuum.absorption_coefficient(temperature, pressure,
                                                                             vmr, grid)*n
            except TypeError:
                k["continua"][mechanism] += continuum.absorption_coefficient(grid)*n
        return k


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
            atmosphere: an xarray Dataset.
            grid: Wavenumber grid array [cm-1].

        Returns:
            An xarray Dataset of absorption coefficients.
        """
        t = atmosphere["temperature"]
        dims = list(t.dims) + ["mechanism", "wavenumber"]
        sizes = tuple([x for x in t.sizes.values()] + [4, grid.size])
        beta = {"{}_absorption".format(name): DataArray(zeros(sizes), dims=dims)
                for name in self.gas.keys()}
        for name, gas in self.gas.items():
            for i in range(t.data.size):
                vmr = {x: atmosphere["vmr_{}".format(x)].data.flat[i]
                       for x in self.gas.keys()}
                k = gas.absorption_coefficient(t.data.flat[i],
                                               atmosphere["pressure"].data.flat[i],
                                               vmr, grid)
                i = unravel_index(i, t.data.shape)
                for j, (source, data) in enumerate(k.items()):
                    indices = tuple(list(i) + [j, slice(None)])
                    if source == "cia":
                        bsum = zeros(grid.size)
                        for values in data.values():
                            bsum += values[:]
                        beta["{}_absorption".format(name)].values[indices] = bsum[:]
                    elif source == "continua":
                        bsum = zeros(grid.size)
                        for values in data.values():
                            bsum += values[:]
                        beta["{}_absorption".format(name)].values[indices] = bsum[:]
                    else:
                        beta["{}_absorption".format(name)].values[indices] = data[:]
        return Dataset(beta)


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
