from time import time

from netCDF4 import Dataset

from PyLBL import Atmosphere, Spectroscopy

from default_atmosphere import grid, pressure, temperature, volume_mixing_ratio


def main():
    timings = {}
    molecules = volume_mixing_ratio.keys()

    tstart = time()
    spectroscopy = Spectroscopy(molecules, "test.db")
    timings["create_database"] = time() - tstart

    tstart = time()
    spectroscopy.load_spectral_inputs()
    timings["read_database"] = time() - tstart

    atmos = Atmosphere(temperature, pressure, volume_mixing_ratio)
    tstart = time()
    absorption_coefficient = spectroscopy.compute_absorption(atmos, grid)
    timings["compute_absorption"] = time() - tstart
    with Dataset("results.nc", "w") as dataset:
        for name, units, data in zip(["pressure", "wavenumber"], ["Pa", "cm-1"],
                                     [atmos.pressure, grid]):
            dataset.createDimension(name, data.size)
            v = dataset.createVariable(name, "f8", (name,))
            v.setncattr("units", units)
            v[:] = data[:]
        for name, data in absorption_coefficient.items():
            v = dataset.createVariable("{}_absorption_coefficient".format(name), "f8",
                                       ("pressure", "wavenumber"))
            v.setncattr("units", "m-1")
            v[...] = data[...]

    print("Timing results:")
    for key, value in timings.items():
        print("{:32s}  -  {}".format(key, value))


if __name__ == "__main__":
    main()
