from logging import getLogger
from re import match
from sqlite3 import connect
from urllib.request import urlopen

from netCDF4 import Dataset
from numpy import asarray, copy, float32, float64, transpose, searchsorted


info = getLogger(__name__).info
SQL_TYPES = {float64: "REAL", int: "INTEGER", str: "TEXT"}
TIPS_REFERENCE_TEMPERATURE = 296.


def ascii_table_records(response, block_size=512):
    """Reads the next line from an ascii table.

    Args:
        response: A http.client.HTTPResponse object.
        block_size: Size in bytes of blocks that will be read.

    Yields:
        Record of the HITRAN database.
    """
    record = ""
    while True:
        block = response.read(block_size).decode("utf-8")
        lines = block.split("\n")
        if lines[-1] == "":
            # If a block ends with a new line character, delete the last
            # element of the list because it will be an empty string.
            del lines[-1]
        for line in lines[:-1]:
            # Return complete lines within a block.
            record += line
            yield record
            record = ""
        if len(block) != block_size:
            # This is the last block.
            yield record + lines[-1]
            break
        elif block.endswith("\n"):
            # No carry-over data between blocks.
            yield lines[-1]
            record = ""
        else:
            # Carry partial last line over to next block.
            record = lines[-1]


def scrub(string):
    """Scrubs user-provided string to prevent database injection.

    Args:
        string: User provided database input.

    Returns:
        The string scrubbed of any trailing spaces, punctuation, or additional text.
    """
    return match(r"([A-Za-z0-9+_-]+)", string.strip()).group(1)


class MoleculeNotFound(Exception):
    pass


class TotalPartitionFunction(object):
    """Total partition function, using data from TIPS 2017 (doi: 10.1016/j.jqsrt.2017.03.045).

    Attributes:
        data: Numpy array of total partition function values (isotopologue, temperature).
        molecule: Molecule chemical formula.
        temperature: Numpy array of temperatures [K].
    """
    def __init__(self, molecule, database=None):
        self.molecule = molecule
        if database is None:
            self.download_from_web()
        else:
            self.load_from_database(database)

    def create_database(self, database):
        """Creates/ingests data into a SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            name = "{}_tips".format(scrub(self.molecule))
            data = transpose(self.data)
            columns = ", ".join(["temperature REAL"] +
                                ["Q_{} REAL".format(i+1) for i in range(data.shape[1])])
            cursor.execute("CREATE TABLE {} ({})".format(name, columns))
            value_subst = ", ".join(["?" for _ in range(data.shape[1] + 1)])
            for i, t in enumerate(self.temperature):
                values = tuple(float64(x) for x in ([t] + data[i, :].tolist()))
                cursor.execute("INSERT INTO {} VALUES ({})".format(name, value_subst), values)
            connection.commit()

    def download_from_web(self):
        """Downloads the data from the internet."""
        url = "http://faculty.uml.edu/Robert_Gamache/Software/temp/Supplementary_file.txt"
        info("Downloading TIPS 2017 data for {} from {}.".format(self.molecule, url))
        self.parse_records(self.records(urlopen(url), self.molecule))

    @property
    def isotopologue(self):
        return [x for x in range(self.data.shape[0])]

    def load_from_database(self, database):
        """Loads data from a previously created SQLite database.

        Args:
            database: Path to SQLite database.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            name = "{}_tips".format(scrub(self.molecule))
            cursor.execute("SELECT * from {}".format(name))
            self.parse_records(cursor.fetchall())

    def parse_records(self, records):
        """Parses all database records and stores the data.

        Args:
            records: A list of iterable database record values.
        """
        temperature, q = [], []
        for record in records:
            temperature.append(record[0])
            q.append(record[1:])
        self.temperature = asarray(temperature, dtype=float32)
        self.data = transpose(asarray(q, dtype=float32))
        info("Found data for {} isotopologues at {} temperatures.".format(*self.data.shape))

    def read_from_dataset(self, path):
        """Reads in total partition function table.

        Args:
            path: Path to input file.
        """
        info("Reading TIPS 2017 data from dataset {}.".format(path))
        with Dataset(path, "r") as dataset:
            self.temperature = copy(dataset.variables["temperature"])
            self.data = copy(dataset.variables["total_partition_function"])

    @staticmethod
    def records(response, molecule):
        """Parses the HTTP table for all records related to the input molecule.

        Args:
            response: A http.client.HTTPResponse object.
            molecule: Molecule id.

        Yields:
            A list of floats from a record from the http table.

        Raises:
            MoleculeNotFound: Failed to find the input molecule.
        """
        found_molecule = False
        num_isotopologues = 0
        for line in ascii_table_records(response):
            if found_molecule:
                if match(r"\s*[A-Za-z0-9]+$", line):
                    break
                elif num_isotopologues > 0:
                    yield [float32(x.strip()) for x in line.split()[:(num_isotopologues+1)]]
                elif match(r"\s*T / K", line):
                    num_isotopologues = sum(x == "Q" for x in line)
            elif line.startswith("c"):
                # Ignore comments.
                continue
            else:
                found_molecule = match(r"\s*{}$".format(molecule), line)
        if not found_molecule:
            raise MoleculeNotFound("molecule {} not found in TIPS 2017 tables.".format(molecule))

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

    def write_to_netcdf(self, path):
        """Writes data to a netCDF dataset.

        Args:
            path: Name of the netCDF4 file that will be created.
        """
        info("Writing TIPS 2017 data to dataset {}.".format(path))
        with Dataset(path, "w") as dataset:
            dataset.createDimension("temperature", self.temperature.size)
            dataset.createDimension("isotopologue", self.data.shape[0])
            v = dataset.createVariable("temperature", float32, dimensions=("temperature",))
            v.setncattr("units", "K")
            v[:] = self.temperature[:]
            v = dataset.createVariable("total_partition_function", float32,
                                       dimensions=("isotopologue", "temperature"))
            v[:, :] = self.data[:, :]
