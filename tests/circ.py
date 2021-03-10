from netCDF4 import Dataset
from numpy import asarray, log, zeros


air_molar_mass = 0.0289647  # Molar mass of air [kg mol-1].
Na = 6.02214076e23  # Avogadro's number [mol-1].
kb = 1.38064852e-23  # Boltzmann constant [J K-1].
g = 9.80  # Acceleration due to gravity [m s-2].


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


class Circ(object):
    def __init__(self, path):
        mb_to_pa = 100.
        with Dataset(path, "r") as dataset:
            self.pressure = dataset.variables["level_pressure"][:]*mb_to_pa
            self.temperature = dataset.variables["level_temperature"][:]
            self.layer_temperature = dataset.variables["layer_temperature"][:]
            self.surface_temperature = dataset.variables["surface_temperature"][:]
            self.vmr = {}
            for molecule in ["CH4", "CO", "CO2", "H2O", "N2O", "O2", "O3"]:
                self.vmr[molecule] = zeros(self.pressure.size)
                self.vmr[molecule][1:] = dataset.variables["{}_abundance".format(molecule)][:]
                self.vmr[molecule][0] = self.vmr[molecule][1]
        self.vmr["N2"] = asarray(self.pressure.size*[0.78,])
