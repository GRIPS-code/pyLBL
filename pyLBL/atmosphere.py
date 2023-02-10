"""Define how atmospheric inputs are handled."""

from re import match


# Map of molecule standard names to chemical formulae.
_standard_name_to_formula = {
    "carbon_dioxide": "CO2",
    "carbon_monoxide": "CO",
    "cfc11": "CFC11",
    "cfc12": "CFC12",
    "methane": "CH4",
    "nitrogen": "N2",
    "nitrous_oxide": "N2O",
    "oxygen": "O2",
    "ozone": "O3",
    "water_vapor": "H2O",
}


class Atmosphere(object):
    """Atmospheric data container with basic data discover methods.

    Attributes:
        dataset: Input xarray Dataset.
        pressure: xarray DataArray object for pressure [Pa].
        temperature: xarray DataArray object for temperature [K].
        gases: Dictionary of xarray DataArray objects for gas mole fractions [mol mol-1].
    """
    def __init__(self, dataset, mapping=None):
        """Initializes an atmosphere object by reading data from an input xarray Dataset.

        Args:
            dataset: xarray Dataset.
            mapping: User-specified mapping of dataset variable names.
        """
        self.dataset = dataset

        # Find the pressure, temperature and gax mixing ratio variables.
        if mapping is None:
            self.pressure = _find_variable(dataset, "air_pressure")
            self.temperature = _find_variable(dataset, "air_temperature")
            self.gases = {x: y for x, y in _gases(dataset)}
        else:
            self.pressure = dataset[mapping["play"]]
            self.temperature = dataset[mapping["tlay"]]
            self.gases = {x: dataset[y] for x, y in mapping["mole_fraction"].items()}


def _find_variable(dataset, standard_name):
    """Finds a variable in a dataset by its standard name attribute.

    Args:
        dataset: xarray Dataset.
        standard_name: String standard name.

    Returns:
        xarray DataArray object.

    Raises:
        ValueError if standard name is not found in the dataset.
    """
    for var in dataset.data_vars.values():
        try:
            if var.attrs["standard_name"] == standard_name:
                return var
        except KeyError:
            continue
    raise ValueError(f"{standard_name} standard name not found in dataset.")


def _gases(dataset):
    """Finds variables that represent gas mole fractions.

    Args:
        dataset: xarray Dataset.

    Yields:
        Gas name (i.e. "H2O") and xarray DataArray.
    """
    for var in dataset.data_vars.values():
        try:
            m = match("mole_fraction_of_([A-Za-z0-9_]+)?_in_air", var.attrs["standard_name"])
        except KeyError:
            continue
        if m:
            yield _standard_name_to_formula[m.group(1)], var
