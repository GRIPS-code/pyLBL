from os import remove
from os.path import join
from urllib.request import urlopen
from zipfile import ZipFile


url = "https://attachment.rrz.uni-hamburg.de/df514eed/coefficients.zip"


def download(directory, name="tmp.zip"):
    zipped = join(directory, name)
    with urlopen(url) as result, open(zipped, "wb") as dl:
        dl.write(result.read())
    with ZipFile(zipped, "r") as f:
        f.extractall(directory)
    remove(zipped)
