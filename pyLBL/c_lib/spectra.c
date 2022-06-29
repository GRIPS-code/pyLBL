#include <math.h>
#include <stdio.h>

#include "spectral_database.h"
#include "voigt.h"


void spectra(double temperature, double pressure, double abundance,
             LineParameter_t parameter, Tips_t tips, double * v, int n,
             int n_per_v, double * k, int cut_off, int remove_pedestal)
{
    double const vlight = 2.99792458e8; /*Speed of light [m s-1].*/
    double const pa_to_atm = 9.86923e-6; /*[atm pa-1].*/
    double const r2 = 2*log(2)*8314.472;
    double const c2 = 1.4387752;

    double p = pressure*pa_to_atm; /*[atm].*/
    double partial_pressure = p*abundance; /*[atm].*/
    double tfact = 296./temperature;

    /*Pressure shift (often 0).*/
    double nu = parameter.nu + p*parameter.delta_air;

    /*Convert for line width in cm-1 at 296 K and 1 atm.*/
    double gamma = (parameter.gamma_air*(p - partial_pressure) +
                   parameter.gamma_self*partial_pressure)*pow(tfact, parameter.n_air);

    /*Calculate Doppler half-width at half-max HWHM in cm-1.*/
    double alpha = (parameter.nu/vlight)*sqrt(r2*temperature/parameter.mass);

    /*Convert for line strength in cm-1.(mol.cm-2)-1 at 296K.*/
    /*Boltzman factor for lower state energy.*/
    double sb = exp(parameter.elower*c2*(temperature - 296.)/(temperature*296.));

    /*Stimulated emission.*/
    double g = exp((-c2*parameter.nu)/temperature);
    double gref = exp((-c2*parameter.nu)/296.);
    double se = (1. - g)/(1. - gref);

    /*Nonlte calculation of absorption coefficient modifiers.*/
    double sq = total_partition_function(tips, 296., parameter.local_iso_id - 1)/
                total_partition_function(tips, temperature, parameter.local_iso_id - 1);

    /*Line strength.*/
    double sw = parameter.sw*sb*se*sq*0.01*0.01;

    /*Find starting and ending indices for the transition.*/
    int s = (floor(nu) - cut_off - v[0])*n_per_v;
    if (s >= n)
    {
        /*Transition does not need to be calculated.*/
        return;
    }
    if (s < 0)
    {
        s = 0;
    }
    int e = (floor(nu) + cut_off + 1 - v[0])*n_per_v;
    if (e >= n)
    {
        e = n - 1;
    }

    /*Calculate absorption coefficient*/
    voigt(v, s, e, nu, alpha, gamma, sw, k);
    if (remove_pedestal != 0)
    {
        double pedestal = k[s];
        if (k[e] < k[s])
        {
            pedestal = k[e];
        }
        int i;
        for (i=s; i<=e; ++i)
        {
            k[i] -= pedestal;
        }
    }
    return;
}
