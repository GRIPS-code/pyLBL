from numpy import arange, exp, log, power, zeros

from .utils import air_number_density, BandedContinuum, Continuum, dry_air_number_density, \
                   LOSCHMIDT, P0, radiation_term, Spectrum, T0, T273


class OxygenContinuum(BandedContinuum):
    def __init__(self):
        self.bands = [OxygenCIAFundamentalContinuum(self.path),
                      OxygenCIANIRContinuum(self.path),
                      OxygenCIANIR2Continuum(self.path),
                      OxygenCIANIR3Continuum(self.path),
                      OxygenVisibleContinuum(self.path),
                      OxygenHerzbergContinuum(self.path),
                      OxygenUVContinuum(self.path)]


class OxygenCIAFundamentalContinuum(Continuum):
    def __init__(self, path):
        self.data = [Spectrum(path, "o2_f"), Spectrum(path, "o2_t")]

    def spectra(self, temperature, pressure, vmr):
        no2 = dry_air_number_density(pressure, temperature, vmr)*vmr["O2"]
        tau_factor = no2*1.e-20*(pressure/P0)*(T273/temperature)
        rad = radiation_term(self.grid()[:], temperature)
        xktfac = (1./T0) - (1./temperature)  # [K-1].
        factor = (1.e20/LOSCHMIDT)  # [cm3].
        return \
            tau_factor*rad[:]*factor*self.data[0].data[:] * \
            exp(self.data[1].data[:]*xktfac)/self.grid()[:]

    def grid(self):
        return self.data[0].wavenumbers()


class OxygenCIANIRContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_inf1")

    def spectra(self, temperature, pressure, vmr):
        no2 = dry_air_number_density(pressure, temperature, vmr)*vmr["O2"]
        ao2 = 1./0.446
        an2 = 0.3/0.446
        tau_factor = \
            (no2/LOSCHMIDT)*(pressure/P0)*(T273/temperature) * \
            (ao2*vmr["O2"] + an2*vmr["N2"] + vmr["H2O"])
        rad = radiation_term(self.grid()[:], temperature)  # [cm-1].
        return tau_factor*rad[:]*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OxygenCIANIR2Continuum(Continuum):
    def __init__(self, path=None):
        self._grid = arange(9100., 11002., 2.)
        self.data = zeros(self._grid.size)
        hw1 = 58.96
        hw2 = 45.04
        for i in range(self._grid.size):
            dv1 = self._grid[i] - 9375.
            dv2 = self._grid[i] - 9439.
            damp1 = exp(dv1/176.1) if dv1 < 0. else 1.
            damp2 = exp(dv2/176.1) if dv2 < 0. else 1.
            o2inf = 0.31831*(((1.166e-04*damp1/hw1)/(1. + (dv1/hw1)*(dv1/hw1))) +
                             ((3.086e-05*damp2/hw2)/(1. + (dv2/hw2)*(dv2/hw2))))*1.054
            self.data[i] = o2inf/self._grid[i]

    def spectra(self, temperature, pressure, vmr):
        no2 = dry_air_number_density(pressure, temperature, vmr)*vmr["O2"]
        n = air_number_density(pressure, temperature, vmr)
        adjwo2 = (no2/n)*(1./vmr["O2"])*no2*1.e-20*(pressure/P0)*(T0/temperature)
        rad = radiation_term(self.grid()[:], temperature)  # [cm-1].
        return adjwo2*rad[:]*self.data[:]

    def grid(self):
        return self._grid[:]


class OxygenCIANIR3Continuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_inf3")

    def spectra(self, temperature, pressure, vmr):
        no2 = dry_air_number_density(pressure, temperature, vmr)*vmr["O2"]
        tau_factor = (no2/LOSCHMIDT)*(pressure/P0)*(T273/temperature)  # [cm3].
        rad = radiation_term(self.grid()[:], temperature)  # [cm-1].
        return tau_factor*rad[:]*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OxygenVisibleContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_invis")

    def spectra(self, temperature, pressure, vmr):
        no2 = dry_air_number_density(pressure, temperature, vmr)*vmr["O2"]  # [cm-3].
        n = air_number_density(pressure, temperature, vmr)
        adjwo2 = (no2/n)*no2*1.e-20*(pressure/P0)*(T273/temperature)  # [cm-3].
        rad = radiation_term(self.grid()[:], temperature)  # [cm-1].
        factor = 1./(LOSCHMIDT*1.e-20*(55.*T273/T0)*(55.*T273/T0)*89.5)  # [cm3].
        return adjwo2*rad[:]*factor*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OxygenHerzbergContinuum(Continuum):
    def __init__(self, path=None):
        self._grid = arange(36000., 100010., 10.)
        self.data = zeros(self._grid.size)
        for i in range(self._grid.size):
            if self._grid[i] <= 36000.:
                self.data[i] = 0.
            else:
                corr = ((40000. - self._grid[i])/4000.)*7.917e-7 \
                       if self._grid[i] <= 40000. else 0.
                yratio = self._grid[i]/48811.0
                self.data[i] = 6.884e-4*yratio*exp(-69.738*power(log(yratio), 2)) - corr

    def spectra(self, temperature, pressure, vmr):
        no2 = dry_air_number_density(pressure, temperature, vmr)*vmr["O2"]  # [cm-3].
        rad = radiation_term(self.grid()[:], temperature)  # [cm-1].
        factor = 1. + 0.83*(pressure/P0)*(T273/temperature)
        return 1.e-20*no2*rad[:]*factor*self.data[:]/self.grid()[:]

    def grid(self):
        return self._grid[:]


class OxygenUVContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_infuv")

    def spectra(self, temperature, pressure, vmr):
        no2 = dry_air_number_density(pressure, temperature, vmr)*vmr["O2"]  # [cm-3].
        rad = radiation_term(self.grid()[:], temperature)  # [cm-1].
        return 1.e-20*no2*rad[:]*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()
