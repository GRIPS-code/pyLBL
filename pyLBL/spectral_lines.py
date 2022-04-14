from math import floor

from numpy import abs, arange, asarray, concatenate, exp, flip, interp, log, min, pi, power, \
                  searchsorted, sqrt, zeros
from scipy.special import wofz

from .tips import TIPS_REFERENCE_TEMPERATURE, TotalPartitionFunction


avogadro = 6.0221409e23 # Avogadro's constant [mol-1].
boltzmann = 1.380658e-16  # Boltzmann constant [erg K-1].
cm_to_m = 0.01  # [m cm-1].
c2 = -1.4387768795689562  # (hc/k) [K cm].
light_speed = 2.99792458e10  # Speed of light [cm s-1].
pa_to_atm = 9.86923e-6  # [Pa-1 atm].
reference_temperature = 296.  # [K].
sqrt_log2 = sqrt(log(2))
sqrt_2 = sqrt(2)
sqrt_2log2 = sqrt(2*log(2))
sqrt_2pi = sqrt(2*pi)


class SpectralLines(object):
    """HITRAN spectral line parameters.

    Attributes:
        d_air: Transition pressure shift factor.
        en: Lower state energy level.
        gamma_air: Air-broadened halfwidth.
        gamma_self: Self-broadened halfwidth.
        iso: Isotopologue id.
        n_air: Air-broadening temperature dependence power.
        s: Transition strength.
        v: Transition wavenumber.
    """
    parameters = ["d_air", "en", "gamma_air", "gamma_self", "iso", "n_air", "s", "v"]
    def __init__(self, transitions):
        """Converts from HAPI transitions objects to numpy arrays.

        Args:
            transitions: List of Transition objects
        """
        for x in self.parameters:
            setattr(self, x, [])
        for transition in transitions:
            self.d_air.append(transition.delta_air)
            self.en.append(transition.elower)
            self.gamma_air.append(transition.gamma_air)
            self.gamma_self.append(transition.gamma_self)
            self.iso.append(transition.local_iso_id)
            self.n_air.append(transition.n_air)
            self.s.append(transition.sw)
            self.v.append(transition.nu)
        for x in self.parameters:
            setattr(self, x, asarray(getattr(self, x)))


def initial_strength_correction(s, q, t, iso, en, v):
    """Initial correction for the line strengths.

    Args:
        s: Line strength [cm].
        q: TotalPartitionFunction object.
        t: Temperature [K].
        iso: HITRAN isotopologue id.
        en: Transition lower state energy [cm-1].
        v: Transition wavenumber [cm-1].

    Returns:
        Corrected line strengths [cm].
    """
    return s*temperature_correct_line_strength(q, t, iso, en, v)


def brute_force(wavenumber, center, strength, gamma, alpha, remove_pedestal, cut_off=25.):
    """Calculates absorption coefficients [cm2].

    Args:
        wavenumber: Numpy array of wavenumbers [cm-1].
        center: Numpy array of transition wavenumbers [cm-1].
        strength: Numpy array of transition strengths [cm].
        gamma: Numpy array of pressure-broadened halfwidths [cm-1].
        alpha: Numpy array of doppler-broadened halfwidths [cm-1].
        remove_pedestal: Flag for removing the MT-CKD continnum pedestal.
        cut_off: Spectral cut-off from line center [cm-1].

    Returns:
        Numpy array of absorption coefficients [cm2].
    """
    left = searchsorted(wavenumber, center - cut_off, side="left")
    right = searchsorted(wavenumber, center + cut_off, side="right")
    k = zeros(wavenumber.size)
    for i in range(center.size):
        dv = abs(wavenumber[left[i]:right[i]] - center[i])
        k[left[i]:right[i]] += strength[i]*voigt_profile(dv, gamma[i], alpha[i])
        if remove_pedestal:
            k[left[i]:right[i]] -= min(k[left[i]], k[right[i] - 1])
    return k


def line_profile(v, s, q, iso, en, d_air, n_air, n_self, gamma_air, gamma_self,
                 mass, t, p, vmr, wavenumber, cut_off=25., remove_pedestal=False):
    """Calculates absorption coefficients [m2].

    Args:
        v: Numpy array of transition wavenumbers [cm-1].
        s: Numpy array of transition strengths [cm].
        q: TotalPartitionFunction object.
        iso: Numpy array of HITRAN isotopologue ids.
        en: Numpy array of transition lower state energy levels [cm-1].
        d_air: Numpy array of air-broadened transition pressure shifts [atm-1 cm-1].
        n_air: Numpy array of air-broadened temperature dependence powers.
        n_self: Numpy array of self-broadened temperature dependence powers.
        gamma_air: Numpy array of air-broadened halfwidths [cm-1 atm-1].
        gamma_self: Numpy array of self-broadened halfwidths [cm-1 atm-1].
        mass: Molecular mass [g mol-1].
        t: Temperature [K].
        p: Pressure [atm].
        vmr: Volume mixing ratio [mol mol-1].
        wavenumber: Numpy array of wavenumbers [cm-1].
        cut_off: Spectral cut-off from line center [cm-1].
        remove_pedestal: Flag for removing the MT-CKD continnum pedestal.

    Returns:
        Numpy array of absorption coefficients [m2].
    """
    center = pressure_shift_transition_wavenumber(v, d_air, p)
    strength = s / temperature_correct_line_strength(q, t, iso, en, v)
    gamma = pressure_broadened_halfwidth(p, p*vmr, t, n_air, n_self, gamma_air, gamma_self)
    alpha = doppler_broadened_halfwidth(t, mass[iso - 1], v)
    k = brute_force(wavenumber, center, strength, gamma, alpha, remove_pedestal, cut_off)
    return k*cm_to_m*cm_to_m


class Gas(object):
    def __init__(self, formula, mass, transitions):
        """Initializes object.

        Args:
            formula: String chemical formula.
            mass: List of isotopologue masses.
            transitions: List of TransitionTable objects.
        """
        self.mass = asarray(mass)
        self.transitions = SpectralLines(transitions)
        self.q = TotalPartitionFunction(formula)
        self.s = initial_strength_correction(self.transitions.s, self.q,
            TIPS_REFERENCE_TEMPERATURE, self.transitions.iso, self.transitions.en,
            self.transitions.v)

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid,
                               continuum="mt-ckd"):
        """Calculates absorption coefficients for the gas.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid [cm-1].

        Returns:
            Absorption coefficients [m2].
        """
        return line_profile(self.transitions.v, self.s, self.q, self.transitions.iso,
                            self.transitions.en, self.transitions.d_air, 
                            self.transitions.n_air, self.transitions.n_air,
                            self.transitions.gamma_air, self.transitions.gamma_self,
                            self.mass, temperature, pressure*pa_to_atm,
                            volume_mixing_ratio, grid, remove_pedestal=continuum == "mt-cdk")


def pressure_shift_transition_wavenumber(line_center, delta_air, pressure):
    """Pressure-shifts transition wavenumbers.

    Args:
        line_center: Transition wavenumber [cm-1].
        delta_air: Air-broadened transition pressure shifts [atm-1 cm-1].
        pressure: Pressure [atm].

    Returns:
        Pressure-shifted transition wavenumber [cm-1].
    """
    return line_center + delta_air*pressure


def temperature_correct_line_strength(partition_function, temperature, isotopologue_id,
                                      lower_state_energy, line_center):
    """Temperature correction factor for a line strength.

    Args:
        partition_function: TotalPartitionFunction object.
        temperature: Temperature [K].
        isotopologue_id: HITRAN isotopologue id.
        lower_state_energy: Transition lower state energy [cm-1].
        line_center: Transition wavenumber [cm-1].

    Returns:
        Temperature correction factor.
    """
    # Divide-by-zeros may occur for transition wavenumbers close to zero, like those
    # for the O16-O17 isotopologue of O2.
    x = partition_function.total_partition_function(temperature, isotopologue_id)
    y = exp(c2*lower_state_energy/temperature)*(1. - exp(c2*line_center/temperature))
    return x/y


def voigt_profile(dv, pressure_halfwidth, doppler_halfwidth):
    """Calculates a Voigt line profile.

    Args:
        dv: Wavenumber distance from line center [cm-1].
        pressure_halfwidth: Pressure-broadened line half-width [cm-1].
        doppler_halfwidth: Doppler-broadened line half-width [cm-1].

    Returns:
        Voigt line profile broadening [cm].
    """
    sigma = doppler_halfwidth / sqrt_2log2
    return wofz((dv + 1j*pressure_halfwidth) / sigma / sqrt_2).real / sigma / sqrt_2pi


def pressure_broadened_halfwidth(pressure, partial_pressure, temperature,
                                 n_air, n_self, gamma_air, gamma_self):
    """Calculates pressure-broadened line halfwidth.

    Args:
        pressure: Pressure [atm].
        partial_pressure: Partial pressure [atm].
        temperature: Temperature [K].
        n_air: Air-broadened temperature dependence power.
        n_self: Self-broadened temperature dependence power.
        gamma_air: Air-broadened halfwidth [cm-1 atm-1].
        gamma_self: Self-broadened halfwidth [cm-1 atm-1].

    Returns:
        Pressure-broadened line halfwidth [cm-1].
    """
    return power((296./temperature), n_air)*(gamma_air*(pressure - partial_pressure)) + \
           power((296./temperature), n_self)*gamma_self*partial_pressure


def doppler_broadened_halfwidth(temperature, mass, transition_wavenumber):
    """Calculate the doppler-broadened line halfwidth.

    Args:
        temperature: Temperature [K].
        mass: Molecular mass [g].
        transition_wavenumber: Transition wavenumber [cm-1].

    Returns:
        Doppler-broadened line halfwidth [cm-1].
    """
    m = mass/avogadro  # Mass/Avogadro's number [g].
    return sqrt_log2*transition_wavenumber*sqrt(2.*boltzmann*temperature/ \
        (m*light_speed*light_speed))
