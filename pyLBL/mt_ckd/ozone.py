from .utils import BandedContinuum, Continuum, dry_air_number_density, radiation_term, \
                   Spectrum, T273


class OzoneContinuum(BandedContinuum):
    def __init__(self):
        self.bands = [OzoneChappuisWulfContinuum(self.path),
                      OzoneHartleyHugginsContinuum(self.path),
                      OzoneUVContinuum(self.path)]


class OzoneChappuisWulfContinuum(Continuum):
    """Ozone continuum in the Chappuis and Wulf band.

    Attributes:
        data: List of Spectrum objects containing data read from an input dataset.
    """
    def __init__(self, path):
        self.data = [Spectrum(path, "x_o3"), Spectrum(path, "y_o3"), Spectrum(path, "z_o3")]

    def spectra(self, temperature, pressure, vmr):
        no3 = dry_air_number_density(pressure, temperature, vmr)*vmr["O3"]
        dt = temperature - T273
        rad = radiation_term(self.grid()[:], temperature)
        return 1.e-20*no3*rad*(self.data[0].data[:] + self.data[1].data[:]*dt +
                               self.data[2].data[:]*dt*dt)/self.grid()[:]

    def grid(self):
        return self.data[0].wavenumbers()


class OzoneHartleyHugginsContinuum(Continuum):
    """Ozone Hartly-Huggins continuum cros sections.

    Attributes:
        data: List of Spectrum objects containing data read from an input dataset.
    """
    def __init__(self, path):
        self.data = [Spectrum(path, "o3_hh0"), Spectrum(path, "o3_hh1"),
                     Spectrum(path, "o3_hh2")]

    def spectra(self, temperature, pressure, vmr):
        no3 = dry_air_number_density(pressure, temperature, vmr)*vmr["O3"]
        dt = temperature - T273
        rad = radiation_term(self.grid()[:], temperature)
        return \
            1.e-20*no3*rad*(self.data[0].data[:]/self.grid()[:]) * \
            (1. + self.data[1].data[:]*dt + self.data[2].data[:]*dt*dt)

    def grid(self):
        return self.data[0].wavenumbers()


class OzoneUVContinuum(Continuum):
    """Ozone ultra-violet continuum coefficients.

    Attributes:
        data: A Spectrum object containing data read from an input dataset.
    """
    def __init__(self, path):
        self.data = Spectrum(path, "o3_huv")

    def spectra(self, temperature, pressure, vmr):
        no3 = dry_air_number_density(pressure, temperature, vmr)*vmr["O3"]
        rad = radiation_term(self.grid()[:], temperature)
        return no3*rad*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()
