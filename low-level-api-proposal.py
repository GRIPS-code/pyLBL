# Spectral lines API:
class SpectralLinesAbstractBaseClass(object):
    """Line-by-line gas optics."""
    def __init__(self, molecule, **kwargs)
        """Initialize object.

        Args:
            molecule: String chemical formula (i.e. "H2O").
            kwargs: Extra optional arguments to configure specific models.
        """
        raise NotImplementedError("you must override this.")

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio,
                               grid, **kwargs):
        """Calculates absorption coefficients for the gas using line-by-line method.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid array [cm-1].
            kwargs: Extra optional arguments to configure specific models.

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        raise NotImplementedError("you must override this.")

# Concrete examples:
# pyrad
class Gas(SpectralLinesAbstractBaseClass):
    """pyrad line-by-line gas optics."""
    def __init__(self, molecule, hitran_database=None, isotopologues=None,
                 line_profile=Voigt(), tips_database=None):
        """Initializes object.

        Args:
            molecule: String chemical formula (i.e. "H2O").
            hitran_database: Path to sqlite hitran database.
            isotopologues: Lists of Isotopologue objects.
            line_profile: Doppler, Lorentz, or Voigt object.
            tips_database: Path to sqlite tips database.
        """
        pass

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio,
                               grid, line_cut_off=25.):
        """Calculates absorption coefficients for the gas using line-by-line method.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid array [cm-1].
            line_cut_off: Cut-off from spectral line center [cm-1].

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        pass

# pygrt
class Gas(SpectralLinesAbstractBaseClass):
    """grtcode line-by-line gas optics."""
    def __init__(self, molecule, hitran_database=None, isotopologues=None,
                 line_profile=Voigt(), tips_database=None, device="host"):
        """Initializes object.

        Args:
            molecule: String chemical formula (i.e. "H2O").
            hitran_database: Path to sqlite hitran database.
            isotopologues: Lists of Isotopologue objects.
            line_profile: Doppler, Lorentz, or Voigt object.
            tips_database: Path to sqlite tips database.
            device: Device to run on.
        """
        pass

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid):
        """Calculates absorption coefficients for the gas using GRTCODE.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid array [cm-1].

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        pass

# pyarts
class PyArtsGas(SpectralLinesAbstractBaseClass):
    """arts line-by-line gas optics."""
    def __init__(self, molecule, workspace=None):
        """Initializes object.

        Args:
            molecule: String chemical formula (i.e. "H2O").
            workspace: Workspace object.  If None, a new one is created.
        """
        pass

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid):
        """Calculates absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid array [cm-1].

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        pass


# Absorption cross section API.  We can make this match the spectral lines API and use
# duck typing.
class SpectralCrossSectionAbstractBaseClass(object):
    """Absorption using HITRAN-supplied absorption cross section tables."""
    def __init__(self, molecule, **kwargs)
        """Initialize object.

        Args:
            molecule: String chemical formula (i.e. "H2O").
            kwargs: Extra optional arguments to configure specific models.
        """
        raise NotImplementedError("you must override this.")

    def absorption_coefficient(self, temperature, pressure, grid, **kwargs):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            grid: Wavenumber grid array [cm-1].
            kwargs: Extra optional arguments to configure specific models.

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        raise NotImplementedError("you must override this.")

# Concrete examples:
# pyrad
class HitranCrossSection(SpectralCrossSectionAbstractBaseClass):
    """pyrad frontend for HITRAN absorption cross-section tables."""
    def __init__(self, molecule, database=None):
        """Downloads HITRAN absorption cross section data tables from the web.

        Args:
            molecule: String chemical formula.
            database: Path to database (not supported yet).
        """
        pass

    def absorption_coefficient(self, temperature, pressure, grid):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            grid: Wavenumber grid array [cm-1].

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        pass


# Collision-induced absorption API.
class CollisionInducedAbsorptionAbstractBaseClass(object):
    """Absorption using HITRAN-supplied absorption cross section tables."""
    def __init__(self, molecule, broadener, **kwargs)
        """Initialize object.

        Args:
            molecule: String chemical formula (i.e. "H2O").
            broadener: String chemical formula (i.e. "H2O") for broadener.
            kwargs: Extra optional arguments to configure specific models.
        """
        raise NotImplementedError("you must override this.")

    def absorption_coefficient(self, temperature, grid, **kwargs):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            temperature: Temperature [K].
            grid: Wavenumber grid array [cm-1].
            kwargs: Extra optional arguments to configure specific models.

        Returns:
            Array of absorption coefficients [m4] on input wavenumber grid.
        """
        raise NotImplementedError("you must override this.")

# Concrete examples:
# pyrad
class HitranCIA(CollisionInducedAbsorptionAbstractBaseClass):
    """pyrad frontent for HITRAN collision-inducd absorption tables."""
    def __init__(self, molecule, broadener, database=None):
        """Initializes object.

        Args:
            molecule: String chemical formula (i.e. "H2O").
            broadener: String chemical formula (i.e. "H2O") for broadener.
            database: Path to database (not supported yet).
        """
        pass

    def absorption_coefficient(self, temperature, grid):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            temperature: Temperature [K].
            grid: Spectral grid array [cm-1].

        Returns:
            Array of absorption coefficients [m4] on input wavenumber grid.
        """
        pass


# Water vapor continuum API.
class WaterVaporContinuumAbstractBaseClass(object):
    """Water vapor continnum."""
    def __init__(self, **kwargs)
        """Initialize object.

        Args:
            kwargs: Extra optional arguments to configure specific models.
        """
        raise NotImplementedError("you must override this.")

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio,
                               grid, **kwargs):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid array [cm-1].
            kwargs: Extra optional arguments to configure specific models.

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        raise NotImplementedError("you must override this.")

# Concrete examples:
# pyrad
class WaterVaporSelfContinuum(WaterVaporContinuumAbstractBaseClass):
    """pyrad frontend for calculating the water vapor self-continuum."""
    def __init__(self, database=None)
        """Initialize object.

        Args:
            database: Path to database of continuum coefficients.
        """
        pass

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Wavenumber grid array [cm-1].
            kwargs: Extra optional arguments to configure specific models.

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        pass


# Ozone continuum API.
class OzoneContinuumAbstractBaseClass(object):
    """Ozone continnum."""
    def __init__(self, **kwargs)
        """Initialize object.

        Args:
            kwargs: Extra optional arguments to configure specific models.
        """
        raise NotImplementedError("you must override this.")

    def absorption_coefficient(self, grid, **kwargs):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            grid: Wavenumber grid array [cm-1].
            kwargs: Extra optional arguments to configure specific models.

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        raise NotImplementedError("you must override this.")

# Concrete examples:
# pyrad
class OzoneContinuum(OzoneContinuumAbstractBaseClass):
    """pyrad frontend for calculating the ozone continuum."""
    def __init__(self, database=None)
        """Initialize object.

        Args:
            database: Path to database of continuum coefficients.
        """
        pass

    def absorption_coefficient(self, grid):
        """Calculates absorption coefficients by interpolation/extrapolation.

        Args:
            grid: Wavenumber grid array [cm-1].

        Returns:
            Array of absorption coefficients [m2] on input wavenumber grid.
        """
        pass
