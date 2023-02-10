"""Controls the TIPS (total internal partition function) calculations."""

from numpy import searchsorted


TIPS_REFERENCE_TEMPERATURE = 296.  # TIPS reference temperature [K].


class TotalPartitionFunction(object):
    """Total partition function, using data from TIPS 2017 (doi: 10.1016/j.jqsrt.2017.03.045).

    Attributes:
        data: Numpy array of total partition function values (isotopologue, temperature).
        molecule: String molecule chemical formula.
        temperature: Numpy array of temperatures [K].
    """
    def __init__(self, molecule, temperature, data):
        self.molecule = molecule
        self.temperature = temperature
        self.data = data

    @property
    def isotopologue(self):
        return [x for x in range(self.data.shape[0])]

    def total_partition_function(self, temperature, isotopologue):
        """Interpolates the total partition function values from the TIPS 2017 table.

        Args:
            temperature: Temperature [K].
            isotopologue: Isotopologue id.

        Returns:
            Total partition function.
        """
        i = isotopologue - 1
        j = searchsorted(self.temperature, temperature, side="left") - 1
        return self.data[i, j] + (self.data[i, j+1] - self.data[i, j]) * \
            (temperature - self.temperature[j])/(self.temperature[j+1] - self.temperature[j])
