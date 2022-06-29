#ifndef SPECTRA_H_
#define SPECTRA_H_

#include "spectral_database.h"


void spectra(double temperature, double pressure, double abundance,
             LineParameter_t parameter, Tips_t tip, double * v, int n,
             int n_per_v, double * k, int cut_off, int remove_pedestal);


#endif
