#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sqlite3.h"

#include "spectra.h"
#include "spectral_database.h"


#define check(err) { \
    if (err == 1) { \
        return 1; \
    } \
}


/*Calculate absorption coefficient spectra.*/
int absorption(double pressure, /*Pressure [Pa].*/
               double temperature, /*Temperature [K].*/
               double volume_mixing_ratio, /*Volume mixing ratio [mol mol-1].*/
               int v0, /*Spectral grid lower bound (inclusive) [cm-1].*/
               int vn, /*Spectral grid upper bound (inclusive) [cm-1].*/
               int n_per_v, /*Number of spectral grid points per wavenumber.*/
               double * k, /*Absorption coefficient [m-1].*/
               char * database, /*Path to the database file.*/
               char * formula, /*Molecue chemical formula.*/
               int cut_off, /*Cut off from line center [cm-1].*/
               int remove_pedestal /*Flag for removing the pedestal.*/
              )
{
    /*Spectral grid.*/
    double dv = 1./n_per_v;
    int n = (vn - v0)*n_per_v;
    double * v = malloc(sizeof(double)*n);
    int i;
    for (i=0; i<n; ++i)
    {
        v[i] = v0 + i*dv;
    }
    memset(k, 0, sizeof(double)*n);

    /*Connect to database.*/
    sqlite3 * connection;
    check(open_database(database, &connection));

    /*Get the molecule id.*/
    int id;
    check(molecule_id(connection, formula, &id));

    /*Read tip data.*/
    Tips_t tips;
    int err = tips_data(connection, id, &tips);
    check(err);
    if (err == -1)
    {
        /*No tips data, so lines can't be calculated.*/
        return 0;
    }

    /*Read the isotopologue masses.*/
    double mass[32];
    memset(mass, 0, sizeof(double)*32);
    check(mass_data(connection, id, mass, 32));

    /*Read the HITRAN line parameters.*/
    char query[256];
    snprintf(query, 256,
             "select nu, sw, gamma_air, gamma_self, n_air, elower, delta_air, "
                 "local_iso_id from transition where molecule_id == %d",
             id);
    sqlite3_stmt * statement;
    compile_statement(connection, query, &statement);

    /*Loop through the lines.*/
    while (sqlite3_step(statement) != SQLITE_DONE)
    {
        LineParameter_t parameter;
        line_parameters(statement, &parameter, mass);
        if (parameter.nu > vn + cut_off + 1 || parameter.nu < v0 - (cut_off + 1))
        {
            break;
        }
        spectra(temperature, pressure, volume_mixing_ratio, parameter, tips,
                v, n, n_per_v, k, cut_off, remove_pedestal);
    }

    /*Clean up clear statement so another query can be made.*/
    sqlite3_finalize(statement);

    /*Close connection to the database.*/
    check(close_database(connection));

    /*Clean up.*/
    free(v);
    free(tips.temperature);
    free(tips.data);
    return 0;
}
