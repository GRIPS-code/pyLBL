from json import loads
from urllib.error import HTTPError
from urllib.request import build_opener, install_opener, ProxyHandler, urlopen


class WebApi(object):
    """Controls access to HITRAN's web API.

    Attributes:
        api_key: String hitran.org api key.
        api_version: String version of the api to use.
        cross_section_directory: String directory where cross section files are located.
        host: URL to retrieve the data from.
        parameters: List of Struct objects describing the HITRAN parameters.
        proxy:
        timestamp: String time stamp telling when the server was accessed.
        transition_directory: String directory where transitions files are located.
    """
    def __init__(self, api_key, api_version="v2", host="https://hitran.org",
                 proxy=None):
        """Initializes object.

        Args:
            api_key: String hitran.org api key.
            api_version: String version of the api to use.
            host: URL to retrieve the data from.
            proxy:
        """
        self.api_key = api_key
        self.api_version = api_version
        self.host = host
        self.proxy = proxy
        server_info = self._download_server_info()
        self.transition_directory = server_info["content"]["data"]["results_dir"]
        self.cross_section_directory = server_info["content"]["data"]["xsec_dir"]
        self.timestamp = server_info["timestamp"]
        self.parameters = self._download_parameters_metadata()

    def _download(self, url, chunk):
        """Downloads data from a url.

        Args:
            url: URL to retrieve data from.
            chunk: Size of data reads in bytes.

        Returns:
            String containing the response data.
        """
        if self.proxy:
            install_opener(build_opener(ProxyHandler(self.proxy)))
        response = urlopen(url)
        data = []
        while True:
            buf = response.read(chunk)
            if not buf:
                break
            data.append(buf.decode("utf-8"))
        return "".join(data)

    def _download_file(self, prefix, name, chunk=64*1024*1024):
        """Downloads a data file from hitran.org.

        Returns:
            String containing the contents of the data file.
        """
        return self._download("/".join([self.host, prefix, name]), chunk)

    def _download_parameters_metadata(self, pattern=None):
        """Downloads metadata about the available HITRAN parameters.

        Args:
            pattern: Substring pattern.

        Returns:
            List of Struct objects containing the response data.
        """
        query = None if pattern is None else Query(name__icontains=pattern)
        return [Struct(**x) for x in
                self._download_section("parameter-metas", query)["content"]["data"]]

    def _download_section(self, api_section, query=None, chunk=1024*1024):
        """Downloads data from the hitran.org website.

        Args:
            api_section: String name of the section of the database to use.
            query: String describing the HTML options.
            chunk: Size of data reads in bytes.

        Returns:
            JSON string containing the response data.
        """
        url = "/".join([self.host, "api", self.api_version, self.api_key, api_section])
        if query is not None:
            url = "?".join([url, query.string])
        if self.proxy:
            install_opener(build_opener(ProxyHandler(self.proxy)))
        return loads(self._download(url, chunk))

    def _download_server_info(self):
        """Downloads internal information about the hitran.org server.

        Returns:
            JSON string containing the response data.
        """
        return self._download_section("info")

    def download_data_sources(self, ids=None):
        """Downloads information about the source of the line data (papers, etc.)

        Args:
            ids: Isotopologue ids.

        Returns:
            JSON string containing the response data.
        """
        query = None if ids is None else Query(id__in=ids)
        return self._download_section("sources", query)["content"]["data"]

    def download_molecules(self):
        """Downloads the molecules available in HITRAN.

        Returns:
            List of Struct objects containing the response data.
        """
        return [Struct(**x) for x in self._download_section("molecules")["content"]["data"]]

    def download_isotopologues(self, molecules):
        """Downloads the isotopologues available in HITRAN.

        Args:
            molecules: List of Struct objects.

        Returns:
            List of Struct objects containing the response data.
        """
        if type(molecules) not in [list, tuple]:
            molecules = [molecules, ]
        ids = [x.id for x in molecules]
        return [Struct(**x) for x in self._download_section("isotopologues",
            Query(molecule_id__in=ids))["content"]["data"]]

    def download_transitions(self, isotopologues, numin, numax, parameters=None):
        """Downloads transitions for isotopologues available in HITRAN.

        Args:
            isotopologues: List of Struct objects.
            numin: Wavenumber lower bound [cm-1].
            numax: Wavenumber upper bound [cm-1].
            parameters: List of parameters to download.

        Returns:
            List of Struct objects containing the response data.
        """
        if type(isotopologues) not in [list, tuple]:
            isotopologues = [isotopologues, ]
        ids = [x.id for x in isotopologues]
        if not ids:
            raise NoIsotopologueError("no isotopologues present.")
        if parameters is None:
            parameters = [x.name for x in self.parameters][:22]
        query = Query(iso_ids_list=ids, numin=numin, numax=numax, head=False,
                      fixwidth=0, request_params=",".join(parameters))
        try:
            name = self._download_section("transitions", query)["content"]["data"]
        except HTTPError:
            raise NoTransitionsError("no transitions found for {}.".format(
                isotopologues[0].molecule_alias))
        data = self._download_file(self.transition_directory, name)

        # Parse the file.
        transitions = []
        type_mapping = {"float": float, "int": int, "str": str}
        types = [type_mapping[x.type] for x in self.parameters]
        for line in data.split("\n"):
            line = line.strip()
            if not line:
                continue
            try:
                transitions.append(Struct(**{x: y(z) for x, y, z in
                                   zip(parameters, types, line.split(","))}))
            except ValueError:
                print("skipping transition: {}".format(line))
        return transitions

    def download_cross_sections(self, molecules):
        """Downloads cross-sections for molecules available in the HITRAN database.

        Args:
            molecules: List of Struct objcts.

        Returns:
            List of Struct objects containing the response data.
        """
        if type(molecules) not in [list, tuple]:
            molecules = [molecules, ]
        ids = [x.id for x in molecules]
        query = Query(molecule_id__in=ids)
        bands = self._download_section("cross-sections", query)["content"]["data"]
        cross_sections = []
        for band in bands:
            data = self._download_file(self.cross_section_directory, band["filename"])
            attrs = {"data": data}
            attrs.update(band)
            cross_sections.append(Struct(**attrs))
        return cross_sections


class NoCrossSectionError(BaseException):
    pass


class NoIsotopologueError(BaseException):
    pass


class NoTransitionsError(BaseException):
    pass


class Query(object):
    """URL parameter (query string) helper class.

    Attributes:
        string: String containing URL parameters.
    """
    def __init__(self, **argv):
        self.string = "&".join(["{}={}".format(key, self.process(argv[key])) for key in argv])

    def __and__(self, q):
        q = Query()
        q.string = "&".join([self.string, q.string])
        return q

    @staticmethod
    def process(val):
        """Process argument value and convert to string."""
        if type(val) in [bool, float, int, str]:
            return str(val)
        if type(val) in [list, set, tuple]:
            return ",".join(str(v) for v in val)
        raise TypeError("bad type for query: '{}'".format(val))


class Struct(object):
    def __init__(self, **attrs):
        self.__dict__.update(attrs)
