from os.path import join
from time import time

from netCDF4 import Dataset
from numpy import arange, asarray

from PyLBL import collision_induced_absorption, continua, cross_sections, models, \
                  molecular_lines

from default_atmosphere import grid, pressure, temperature, volume_mixing_ratio


def main():
    timings = {}
    with Dataset("PyLBL-output.nc", "w") as dataset:
        for data, name, units in zip([pressure, grid], ["pressure", "wavenumber"],
                                     ["Pa", "cm-1"]):
            dataset.createDimension(name, data.size)
            v = dataset.createVariable(name, "f8", (name,))
            v.setncattr("units", units)
            v[:] = data[:]

        for model in models:
            # Line-by-line spectra.
            for molecule in ["H2O", "CO2", "O3", "N2O", "CH4", "CO", "O2"]:
                gas = molecular_lines[model](molecule)
                v = dataset.createVariable("{}_{}_lines".format(molecule, model), "f8",
                                           ("pressure", "wavenumber"))
                v.setncattr("units", "m2")
                key = "{}-{}-lines".format(molecule, model)
                timings[key] = 0.
                for layer in range(pressure.size):
                    tstart = time()
                    v[layer, :] = gas.absorption_coefficient(temperature[layer], pressure[layer],
                                                             volume_mixing_ratio[molecule][layer],
                                                             grid)
                    tend = time()
                    timings[key] += tend - tstart

            # Hitran supplied absorption cross sections.
            for molecule in ["CFC-11", "CFC-12"]:
                try:
                    gas = cross_sections[model](molecule)
                except KeyError:
                    continue
                v = dataset.createVariable("{}_{}_xsec".format(molecule, model), "f8",
                                           ("pressure", "wavenumber"))
                v.setncattr("units", "m2")
                key = "{}-{}-xsec".format(molecule, model)
                timings[key] = 0.
                for layer in range(pressure.size):
                    tstart = time()
                    v[layer, :] = gas.absorption_coefficient(temperature[layer], pressure[layer],
                                                             grid)
                    tend = time()
                    timings[key] += tend - tstart

            # Hitran supplied collision-induced absorption.
            for collision in [["O2", "N2"], ["O2", "O2"], ["N2", "N2"]]:
                try:
                    gas = collision_induced_absorption[model](*collision)
                except KeyError:
                    continue
                v = dataset.createVariable("{}_{}_{}_cia".format(*collision, model), "f8",
                                           ("pressure", "wavenumber"))
                v.setncattr("units", "m5")
                key = "{}-{}-{}-cia".format(*collision, model)
                timings[key] = 0.
                for layer in range(pressure.size):
                    tstart = time()
                    v[layer, :] = gas.absorption_coefficient(temperature[layer], grid)
                    tend = time()
                    timings[key] += tend - tstart

            # Continua.
            try:
                absorption_continua = continua[model]
            except KeyError:
                continue
            for molecule, continuum in absorption_continua.items():
                gas = continuum()
                print("continua", molecule, model)
                v = dataset.createVariable("{}_{}_continuum".format(molecule, model),
                                           "f8", ("pressure", "wavenumber"))
                v.setncattr("units", "m2")
                key = "{}-{}-continuum".format(molecule, model)
                timings[key] = 0.
                if molecule.startswith("H2O"): molecule = "H2O"
                for layer in range(pressure.size):
                    tstart = time()
                    try:
                        v[layer, :] = gas.absorption_coefficient(temperature[layer], pressure[layer],
                                                                 volume_mixing_ratio[molecule][layer],
                                                                 grid)
                    except TypeError:
                        v[layer, :] = gas.absorption_coefficient(grid)
                    tend = time()
                    timings[key] += tend - tstart

    print("Timing results:")
    for key, value in timings.items():
        print("{:32s}  -  {}".format(key, value))


if __name__ == "__main__":
    main()
