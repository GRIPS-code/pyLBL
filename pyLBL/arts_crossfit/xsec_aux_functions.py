"""
Created on Fri Sep 25 16:10:11 2020

@author: Manfred Brath

This file contains the functions that are needed for the harmonization of the Hitran
absorption cross section data and for the calculations of the fit coefficients.

Modified by menzel-gfdl.
"""
from numpy import shape, sum, zeros


def calculate_xsec(temperature, pressure, coeffs):
    """
    Low level function to calculate the absorption cross section from the fitted
    coefficients

    Args:
        temperature: Float temperature [K].
        pressure: Float pressure [Pa].
        coeffs: matrix of fit coefficients.

    Returns:
        Vector of absorption cross section [m2].
    """

    # The fit model
    # 2d quadratic fit:
    # z = p00 + p10*x + p01*y + p20*x*x
    # z = xsec
    # x = T
    # y = P

    # coeffs[0,:] = p00
    # coeffs[1,:] = p10
    # coeffs[2,:] = p01
    # coeffs[3,:] = p20

    # distinguish if we calculate xsec for a lot of frequencies
    if len(shape(coeffs)) > 1:
        poly = zeros(4)
        poly[0] = 1
        poly[1] = temperature
        poly[2] = pressure
        poly[3] = temperature*temperature

        # allocate
        xsec = zeros(shape(coeffs))

        for i in range(4):
            xsec[i, :] = coeffs[i, :] * poly[i]

    # or for a lot of states
    else:
        poly = zeros((4, len(temperature)))
        poly[0, :] = 1.
        poly[1, :] = temperature
        poly[2, :] = pressure
        poly[3, :] = temperature*temperature

        # allocate
        xsec = zeros((len(coeffs), len(temperature)))

        for i in range(4):
            xsec[i, :] = coeffs[i] * poly[i, :]

    xsec = sum(xsec, axis=0)

    return xsec


def calculate_xsec_fullmodel(temperature, pressure, coeffs):
    """
    Function to calculate the absorption cross section from the fitted
    coefficients including check for negative values.

    Args:
        temperature: Float temperature [K].
        pressure: Float pressure in [Pa].
        coeffs: matrix fit coefficients.

    Returns:
        Vector of absorption cross section in [m2].
    """

    # The fit model
    # 2d quadratic fit:
    # z = p00 + p10*x + p01*y + p20*x*x

    # z = Xsec
    # x = T
    # y = P

    # coeffs[0,:] = p00
    # coeffs[1,:] = p10
    # coeffs[2,:] = p01
    # coeffs[3,:] = p20

    # calculate raw xsecs
    xsec = calculate_xsec(temperature, pressure, coeffs)

    # Check for negative values and remove them without introducing bias, meaning
    # the integral over the spectrum must not change.
    logic = xsec < 0
    if sum(logic) > 0:

        # original sum over spectrum
        sumx_org = sum(xsec)

        # remove negative values
        xsec[logic] = 0

        if sumx_org >= 0:
            # estimate ratio between altered and original sum of spectrum
            w = sumx_org / sum(xsec)

            # scale altered spectrum
            xsec = xsec * w

    return xsec
