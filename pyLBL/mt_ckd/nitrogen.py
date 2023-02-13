from numpy import power, zeros

from .utils import BandedContinuum, Continuum, dry_air_number_density, LOSCHMIDT, P0, \
                   radiation_term, Spectrum, T0, T273


class NitrogenContinuum(BandedContinuum):
    def __init__(self):
        self.bands = [NitrogenCIAPureRotationContinuum(self.path),
                      NitrogenCIAFundamentalContinuum(self.path),
                      NitrogenCIAFirstOvertoneContinuum(self.path)]


class NitrogenCIAPureRotationContinuum(Continuum):
    def __init__(self, path):
        self.data = {296: [Spectrum(path, "ct_296"), Spectrum(path, "sf_296")],
                     220: [Spectrum(path, "ct_220"), Spectrum(path, "sf_220")]}

    def spectra(self, temperature, pressure, vmr):
        nn2 = dry_air_number_density(pressure, temperature, vmr)*vmr["N2"]
        tau_factor = (nn2/LOSCHMIDT)*(pressure/P0)*(T273/temperature)
        rad = radiation_term(self.grid()[:], temperature)
        factor = (temperature - T0)/(220. - T0)
        c = self.data[296][0].data[:]*power(self.data[220][0].data[:]/self.data[296][0].data[:],
                                            factor)
        s = self.data[296][1].data[:]*power(self.data[220][1].data[:]/self.data[296][1].data[:],
                                            factor)
        fo2 = (s[:] - 1.)*vmr["N2"]/vmr["O2"]
        return tau_factor*rad[:]*c[:]*(vmr["N2"] + fo2[:]*vmr["O2"] + vmr["H2O"])

    def grid(self):
        return self.data[296][0].wavenumbers()


class NitrogenCIAFundamentalContinuum(Continuum):
    def __init__(self, path):
        self.data = [Spectrum(path, "xn2_272"), Spectrum(path, "xn2_228"),
                     Spectrum(path, "a_h2o")]

    def spectra(self, temperature, pressure, vmr):
        nn2 = dry_air_number_density(pressure, temperature, vmr)*vmr["N2"]
        tau_factor = (nn2/LOSCHMIDT)*(pressure/P0)*(T273/temperature)
        rad = radiation_term(self.grid()[:], temperature)

        xtfac = (1./temperature - 1./272.)/(1./228. - 1./272.)
        ao2 = 1.294 - 0.4545*temperature/T0
        c0 = zeros(self.data[0].data.size)
        c0[1: -1] = self.data[0].data[1: -1]*power(self.data[1].data[1: -1] /
                                                   self.data[0].data[1: -1], xtfac)
        c0 = c0[:]/self.grid()[:]
        c1 = ao2*c0[:]
        c2 = (9./7.)*self.data[2].data[:]*c0[:]
        return tau_factor*rad[:]*(c0[:]*vmr["N2"] + vmr["O2"]*c1[:] + vmr["H2O"]*c2[:])

    def grid(self):
        return self.data[0].wavenumbers()


class NitrogenCIAFirstOvertoneContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "xn2")

    def spectra(self, temperature, pressure, vmr):
        nn2 = dry_air_number_density(pressure, temperature, vmr)*vmr["N2"]
        tau_factor = (nn2/LOSCHMIDT)*(pressure/P0)*(T273/temperature) * \
                     (vmr["N2"] + vmr["O2"] + vmr["H2O"])
        rad = radiation_term(self.grid()[:], temperature)
        return tau_factor*rad[:]*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()
