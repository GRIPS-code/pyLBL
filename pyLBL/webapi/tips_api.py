from re import match
from urllib.request import urlopen

from numpy import asarray, float32, transpose


class TipsWebApi(object):
    """Controls access to TIPS web API.

    Attributes:
        url: String url where data is downloaded from.
        timestamp: Dictionary of time stamps from when the molecule data was downloaded.
    """
    def __init__(self):
        """Constructs the object."""
        self.url = "http://faculty.uml.edu/Robert_Gamache/Software/temp/Supplementary_file.txt"
        self.timestamp = {}

    def download(self, molecule):
        """Downloads the data from the internet.

        Args:
            molecule: String molecule chemical formula.

        Returns:
            temperature: Numpy array of temperatures.
            data: Numpy array of data values.
        """
#       self.timestamp[molecule] = 
        return self._parse_records(self._records(urlopen(self.url), molecule))

    @staticmethod
    def _ascii_table_records(response, block_size=512):
        """Reads the next line from an ascii table.

        Args:
            response: A http.client.HTTPResponse object.
            block_size: Size in bytes of blocks that will be read.

        Yields:
            String record of the http table.
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
                try:
                    yield record + lines[-1]
                except IndexError:
                    yield record
                break
            elif block.endswith("\n"):
                # No carry-over data between blocks.
                yield lines[-1]
                record = ""
            else:
                # Carry partial last line over to next block.
                record = lines[-1]

    @staticmethod
    def _parse_records(records):
        """Parses all table records and stores the data.

        Args:
            records: A list of iterable database record values.

        Returns:
            temperature: Numpy array of temperatures.
            data: Numpy array of data values.
        """
        temperature, q = [], []
        for record in records:
            if record:
                temperature.append(record[0])
                q.append(record[1:])
        temperature = asarray(temperature, dtype=float32)
        data = transpose(asarray(q, dtype=float32))
        return temperature, data

    def _records(self, response, molecule):
        """Parses the HTTP table for all records related to the input molecule.

        Args:
            response: A http.client.HTTPResponse object.
            molecule: Molecule id.

        Yields:
            A list of floats from a record from the http table.

        Raises:
            NoMoleculeError: Failed to find the input molecule.
        """
        found_molecule = False
        num_isotopologues = 0
        for line in self._ascii_table_records(response):
            if found_molecule:
                if match(r"\s*[A-Za-z0-9+]+$", line):
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
            raise NoMoleculeError(f"molecule {molecule} not found in TIPS 2017 tables.")


class NoMoleculeError(BaseException):
    """No TIPS data found for this molecule."""
    pass
