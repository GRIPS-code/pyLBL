from numpy import power, zeros

from .utils import air_number_density, BandedContinuum, Continuum, dry_air_number_density, \
                   P0, radiation_term, Spectrum, subgrid_bounds, T0


class WaterVaporSelfContinuum(BandedContinuum):
    def __init__(self):
        self.bands = [WaterVaporARMSelfContinuum(self.path),]


class WaterVaporARMSelfContinuum(Continuum):
    """Water vapor self continuum coefficients.

    Attributes:
        data: Dictionary that maps temperatures (keys) to Spectrum objects containing
              data read from an input dataset (values).
    """
    def __init__(self, path):
        self.data = {296: Spectrum(path, "bs296"),
                     260: Spectrum(path, "bs260")}

    def spectra(self, temperature, pressure, vmr):
        t_factor = (temperature - T0)/(260. - T0)
        nh2o = dry_air_number_density(pressure, temperature, vmr)*vmr["H2O"]
        n = air_number_density(pressure, temperature, vmr)
        rad = radiation_term(self.grid()[:], temperature)
        return \
            nh2o*(nh2o/n)*(pressure/P0)*(T0/temperature)*1.e-20*rad * \
            self.data[296].data[:]*power(self.data[260].data[:]/self.data[296].data[:],
                                         t_factor)

    def grid(self):
        return self.data[296].wavenumbers()


class WaterVaporForeignContinuum(BandedContinuum):
    def __init__(self):
        self.bands = [WaterVaporIASIForeignContinuum(self.path),]


class WaterVaporIASIForeignContinuum(Continuum):
    """Water vapor foreign continuum coefficients.

    Attributes:
        data: Spectrum object containing data read from an input dataset.
        scale: Array of scaling factors.
    """
    def __init__(self, path):
        self.data = Spectrum(path, "bfh2o")
        x = Spectrum(path, "xfac_rhu")
        self.scale = zeros(self.data.data.size)
        lower, upper = subgrid_bounds(self.data.grid, x.grid)
        self.scale[lower + 1:upper + 1] = x.data[1:]
        self.scale[lower] = self.scale[lower + 1]
        u = upper + 1
        w = self.grid()[u:]
        vdelsq1 = (w[:] - 255.67)*(w[:] - 255.67)
        vf1 = power((w[:] - 255.67)/57.83, 8)
        vdelmsq1 = (w[:] + 255.67)*(w[:] + 255.67)
        vmf1 = power((w[:] + 255.67)/57.83, 8)
        vf2 = power(w[:]/630., 8)
        self.scale[u:] = \
            1. + (0.06 - 0.42*((57600./(vdelsq1[:] + 57600. + vf1[:])) +
                  (57600./(vdelmsq1[:] + 57600. + vmf1[:]))))/(1. + 0.3*vf2[:])

    def spectra(self, temperature, pressure, vmr):
        nh2o = dry_air_number_density(pressure, temperature, vmr)*vmr["H2O"]
        n = air_number_density(pressure, temperature, vmr)
        rad = radiation_term(self.grid()[:], temperature)
        return \
            (1. - (nh2o/n))*(pressure/P0)*(T0/temperature)*1.e-20*nh2o*rad * \
            self.scale[:]*self.data.data[:]

    def grid(self):
        return self.data.wavenumbers()
