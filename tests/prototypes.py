from os.path import join
from time import time

from netCDF4 import Dataset
from numpy import arange, asarray

from PyLBL import collision_induced_absorption, continua, cross_sections, models, \
                  molecular_lines


def main():
    #Set up some test atmospheric points.
    grid = arange(1., 3250., 0.01)
    pressure = asarray([117., 1032., 11419., 98388.])  # [Pa].
    temperature = asarray([269.01, 227.74, 203.37, 288.99])  # [K].
    volume_mixing_ratio = {
        "H2O": asarray([5.244536e-06, 4.763972e-06, 3.039952e-06, 6.637074e-03]),
        "CO2": asarray([0.00036, 0.00036, 0.00036, 0.00035999]),
        "O3": asarray([2.936688e-06, 7.415223e-06, 2.609510e-07, 6.859128e-08]),
        "N2O": asarray([1.050928e-08, 1.319584e-07, 2.895416e-07, 3.199949e-07]),
        "CH4": asarray([2.947482e-07, 8.817705e-07, 1.588336e-06, 1.700002e-06]),
        "CO": asarray([3.621464e-08, 1.761450e-08, 3.315927e-08, 1.482969e-07]),
        "O2": asarray([0.209, 0.209, 0.2090003, 0.208996])
    }

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
