from argparse import ArgumentParser
from os import environ

from pyLBL import Database, HitranWebApi


def main(database_path, api_key, molecules=["H2O",]):
    database = Database(database_path)
    webapi = HitranWebApi(api_key)
    database.create(webapi, molecules=molecules)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("path", help="Path to database file.")
    parser.add_argument("--api_key", help="HITRAN api key.", default="")
    args = parser.parse_args()
    if not args.api_key:
        try:
            args.api_key = environ["HITRAN_API_KEY"]
        except KeyError:
            pass
    main(args.path, args.api_key)
