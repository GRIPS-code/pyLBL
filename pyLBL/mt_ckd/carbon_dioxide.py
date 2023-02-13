from numpy import ones, power

from .utils import BandedContinuum, Continuum, dry_air_number_density, P0, radiation_term, \
                   Spectrum, subgrid_bounds, T0


class CarbonDioxideContinuum(BandedContinuum):
    def __init__(self):
        self.bands = [CarbonDioxideHartmannContinuum(self.path),]


class CarbonDioxideHartmannContinuum(Continuum):
    """Carbon dioxide continuum coefficients.

    Attributes:
        data: Spectrum object containing data read from an input dataset.
        t_correction: Array of temperature correction coefficients.
        xfac_co2: Array of chi-factors.
    """
    def __init__(self, path):
        self.data = Spectrum(path, "bfco2")

        x = Spectrum(path, "tdep_bandhead")
        lower, upper = subgrid_bounds(self.data.grid, x.grid)
        self.t_correction = ones(self.data.data.size)
        self.t_correction[lower:upper + 1] = x.data[:]

        x = Spectrum(path, "x_factor_co2")
        lower, upper = subgrid_bounds(self.data.grid, x.grid)
        self.xfac_co2 = ones(self.data.data.size)
        self.xfac_co2[lower:upper + 1] = x.data[:]

    def spectra(self, temperature, pressure, vmr):
        nco2 = dry_air_number_density(pressure, temperature, vmr)*vmr["CO2"]
        rad = radiation_term(self.grid()[:], temperature)
        return \
            nco2*1.e-20*(pressure/P0)*(T0/temperature)*rad[:] * \
            self.xfac_co2[:]*power(temperature/246., self.t_correction[:]) * \
            self.data.data[:]

    def grid(self):
        return self.data.wavenumbers()
