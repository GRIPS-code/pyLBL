from argparse import ArgumentParser

from numpy import arange
from pyLBL import Database, Spectroscopy


def variable(data, units, standard_name):
    return (["z",], data, {"units": units, "standard_name": standard_name})


def atmosphere():
    from xarray import Dataset
    temperature = [288.99,]  # [K].
    pressure = [98388.,]  # [Pa].
    xh2o = [0.006637074,]  # [mol mol-1].
    xcfc11 = [2.783e-10,]  # [mol mol-1].
    return Dataset(
        data_vars={
            "play": variable(pressure, "Pa", "air_pressure"),
            "tlay": variable(temperature, "K", "air_temperature"),
            "xh2o": variable(xh2o, "mol mol-1", "mole_fraction_of_water_vapor_in_air"),
            "xcfc11": variable(xcfc11, "mol mol-1", "mole_fraction_of_cfc11_in_air"),
         },
         coords={
             "layer": (["z",], [1,])
         },
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("path", help="Path to database file.")
    args = parser.parse_args()
    database = Database(args.path)
    dataset = atmosphere()
    grid = arange(1., 5000., 1.)
    s = Spectroscopy(dataset, grid, database)
    print(s.list_molecules())
    output = s.compute_absorption(output_format="all")
    output.to_netcdf("out-all.nc")
