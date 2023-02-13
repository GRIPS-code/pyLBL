from numpy import shape, zeros
from scipy.interpolate import interp1d
from xarray import open_dataset

from .xsec_aux_functions import calculate_xsec_fullmodel


class CrossSection(object):
    def __init__(self, formula, path):
        """Initializes the object.

        Args:
            formula: String chemical formula.
            path: Path to the data file.
        """
        self.formula = formula
        self.path = path

    def absorption_coefficient(self, grid, temperature, pressure):
        """Calculates absorption cross sections.

        Args:
            grid: Numpy array of wavenumbers [cm-1].
            temperature: Temperature [K].
            pressure: Pressure [Pa].

        Returns:
            Numpy array of absorption cross sections in [m2].
        """
        with open_dataset(self.path) as xsec_data:
            # Convert desired wavenumber to frequency [Hz].
            c0 = 299792458.0  # Speed of light [m s-1].
            freq_user = grid * c0 * 100
            xsec_user = zeros(shape(grid))
            bands = xsec_data.bands.data
            for m in bands:
                arg = f"band{m}"
                # frequency of data in [Hz]
                freq_data = xsec_data[arg + "_fgrid"].data.transpose()
                # fit coefficients of band m
                coeffs_m = xsec_data[arg + "_coeffs"].data.transpose()
                # Calculate the cross section on their internal frequency grid
                xsec_temp = calculate_xsec_fullmodel(temperature, pressure, coeffs_m)
                # Interpolate cross sections to user grid
                f_int = interp1d(freq_data, xsec_temp, fill_value=0., bounds_error=False)
                xsec_user_m = f_int(freq_user)
                xsec_user = xsec_user + xsec_user_m
            return xsec_user
