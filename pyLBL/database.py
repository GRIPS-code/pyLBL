"""Manages how the spectral input data is acquired and handled."""

from os import listdir
from os.path import abspath, join
from pathlib import Path
from re import match

from numpy import asarray, reshape
from sqlalchemy import Column, create_engine, Float, ForeignKey, Integer, select, String
from sqlalchemy.orm import declarative_base, Session

from .arts_crossfit.webapi import download
from .tips import TotalPartitionFunction
from .webapi import NoIsotopologueError, NoMoleculeError, NoTransitionsError, TipsWebApi


Base = declarative_base()


def _ingest_molecule_metadata(session, molecule):
    """Stores the molecule metadata in the database.

    Args:
        session: Session object.
        molecule: Struct object.
    """
    session.add(
        MoleculeTable(
            id=molecule.id,
            stoichiometric_formula=molecule.stoichiometric_formula,
            ordinary_formula=molecule.ordinary_formula,
            common_name=molecule.common_name,
        )
    )


def _ingest_molecule_aliases(session, molecule):
    """Stores the molecule aliases in the database.

    Args:
        session: Session object.
        molecule: Struct object.
    """
    aliases = [x["alias"] for x in molecule.aliases]
    for alias in aliases:
        session.add(
            MoleculeAliasTable(
                alias=alias,
                molecule=molecule.id,
            )
        )


def _ingest_isotopologues(session, molecule, hitran_webapi):
    """Stores the isotopologues in the database.

    Args:
        session: Session object.
        molecule: Struct object.
        hitran_webapi: HitranWepApi object.

    Returns:
        List of Struct objects.
    """
    isotopologues = hitran_webapi.download_isotopologues(molecule)
    for isotopologue in isotopologues:
        session.add(
            IsotopologueTable(
                id=isotopologue.id,
                molecule_id=molecule.id,
                isoid=isotopologue.isoid,
                iso_name=isotopologue.iso_name,
                abundance=isotopologue.abundance,
                mass=isotopologue.mass,
            )
        )
    return isotopologues


def _ingest_transitions(session, molecule, isotopologues, hitran_webapi):
    """Stores the transitions in the database.

    Args:
        session: Session object.
        molecule: Struct object.
        isotopologues: List of Struct object.
        hitran_webapi: HitranWepApi object.
    """
    parameters = ["global_iso_id", "molec_id", "local_iso_id", "nu", "sw",
                  "gamma_air", "gamma_self", "n_air", "delta_air", "elower"]
    transitions = hitran_webapi.download_transitions(isotopologues, 0., 1.e8, parameters)
    for transition in transitions:
        session.add(
            TransitionTable(
                global_iso_id=transition.global_iso_id,
                molecule_id=molecule.id,
                local_iso_id=transition.local_iso_id,
                nu=transition.nu,
                sw=transition.sw,
                gamma_air=transition.gamma_air,
                gamma_self=transition.gamma_self,
                n_air=transition.n_air,
                delta_air=transition.delta_air,
                elower=transition.elower,
            )
        )


def _ingest_tips(session, molecule, tips_webapi):
    """Stores the TIPS data in the database.

    Args:
        session: Session object.
        molecule: Struct object.
        tips_webapi: TipsWepApi object.
    """
    temperature, data = tips_webapi.download(molecule.ordinary_formula)
    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            session.add(
                TipsTable(
                    molecule_id=molecule.id,
                    isotopologue_id=x,
                    temperature=temperature[y],
                    data=data[x, y],
                )
            )


class Database(object):
    """Spectral line parameter database.

    Attributes:
        engine: SQLAlchemy Engine object.
    """
    def __init__(self, path, echo=False):
        """Connects to the database and creates tables.

        Args:
            path: Path to database.
            echo: Print database SQL commands.
        """
        self.cross_section_directory = None
        self.engine = create_engine(f"sqlite+pysqlite:///{path}", echo=echo, future=True)
        Base.metadata.create_all(self.engine)
        self.path = path

    def create(self, hitran_webapi, molecules="all", tips_webapi=None,
               cross_section_directory=".cross-sections"):
        """Downloads data from HITRAN and TIPS and inserts it into the database tables.

        Args:
            hitran_webapi: HitranWebApi object.
            molecules: List of string molecule chemical formulae.
            tips_webapi: TipsWebApi object.
            cross_section_directory: Directory to download cross sections to.
        """
        if tips_webapi is None:
            tips_webapi = TipsWebApi()

        with Session(self.engine, future=True) as session:
            all_molecules = hitran_webapi.download_molecules()
            num_molecules = len(all_molecules) if molecules == "all" else len(molecules)
            for i, molecule in enumerate(all_molecules):

                # Support for only using a subset of molecules.
                if molecules != "all" and molecule.ordinary_formula not in molecules:
                    continue

                # Write out progress.
                print("Working on molecule {} / {} ({})".format(
                    i + 1, num_molecules, molecule.ordinary_formula)
                )

                # Store the molecule metadata.
                _ingest_molecule_metadata(session, molecule)

                # Store the molecule aliases.
                _ingest_molecule_aliases(session, molecule)

                # Store isotopologues.
                isotopologues = _ingest_isotopologues(session, molecule, hitran_webapi)

                # Store the transitions.
                try:
                    _ingest_transitions(session, molecule, isotopologues, hitran_webapi)
                except NoIsotopologueError:
                    print(f"No isotopologues for molecule {molecule.ordinary_formula}.")
                    continue
                except NoTransitionsError:
                    print(f"No transitions for molecule {molecule.ordinary_formula}.")
                    continue

                # Store the TIPS data.
                try:
                    _ingest_tips(session, molecule, tips_webapi)
                except NoMoleculeError:
                    print(f"No molecule {molecule.ordinary_formula} found in TIPS database.")
                    continue

                # Commit the updates to the database.
                session.commit()

        # Add in the ARTS Crossfit data files.
        self.cross_section_directory = cross_section_directory
        Path(self.cross_section_directory).mkdir(parents=True, exist_ok=True)
        download(self.cross_section_directory)
        with Session(self.engine, future=True) as session:
            self._ingest_arts_crossfit(session, molecules)
            session.commit()

    def _formula(self, session, molecule_id):
        """Helper function that retrieves a molecule's chemical formula.

        Args:
            session: SQLAlchemy Session object.
            molecule_id: MoleculeTable integer primary key.

        Returns:
            MoleculeTable string chemical formula (i.e., H2O).
        """
        stmt = select(MoleculeTable.ordinary_formula).filter_by(id=molecule_id)
        return session.execute(stmt).all()[0][0]

    def _ingest_arts_crossfit(self, session, molecules):
        """Stores the arts-crossfit data in the database.

        Args:
            session: Session object.
            molecules: List of string molecule chemical formulae.
        """
        dir_path = join(self.cross_section_directory, "coefficients")
        for path in listdir(dir_path):
            # Find chemical formula from file name.
            regex = match(r"([A-Za-z0-9]+).nc", path)
            if not regex:
                continue
            formula = regex.group(1)

            # Support for only using a subset of molecules.
            if molecules != "all" and formula not in molecules:
                continue

            try:
                # Get the id of the molecule in the database.
                molecule_id = self._molecule_id(session, formula)
            except AliasNotFoundError:
                # The molecule doesn't exist in the database, so add it in.
                session.add(
                    MoleculeTable(
                        stoichiometric_formula=formula,
                        ordinary_formula=formula,
                        common_name=formula,
                    )
                )
                session.commit()

                # Query for the molecule id.
                stmt = select(MoleculeTable.id).filter_by(ordinary_formula=formula)
                molecule_id = session.execute(stmt).all()[0][0]
                # Add molecule to the MoleculeAliasTable.
                session.add(
                    MoleculeAliasTable(
                        alias=formula,
                        molecule=molecule_id,
                    )
                )
                session.commit()

            # Store path to the cross section data file.
            full_path = abspath(join(dir_path, path))
            session.add(
                ArtsCrossFitTable(
                    molecule_id=molecule_id,
                    path=full_path,
                )
            )

    def _mass(self, session, molecule_id):
        """Helper function that retrieves a molecule's isotopologue masses.

        Args:
            session: SQLAlchemy Session object.
            molecule_id: MoleculeTable integer primary key.

        Returns:
            List of float isotopologue masses.

        Raises:
            IsotopologuesNotFoundError if no isotopologue is found.
        """
        stmt = select(IsotopologueTable.mass).filter_by(molecule_id=molecule_id)
        result = [x[0] for x in session.execute(stmt).all()]
        if not result:
            raise IsotopologuesNotFoundError(
                      f"isotopologues not found for molecule {molecule_id}."
                  )
        return result

    def _molecule_id(self, session, name):
        """Helper function that retrieves a molecules id.

        Args:
            session: SQLAlchemy Session object.
            name: String molecule alias.

        Returns:
            MoleculeTable integer primary key.

        Raises:
            AliasNotFoundError if no molecule alias is found.
        """
        stmt = select(MoleculeAliasTable.molecule).filter_by(alias=name)
        try:
            return session.execute(stmt).all()[0][0]
        except IndexError:
            raise AliasNotFoundError(f"{name} not found in database.")

    def _transitions(self, session, molecule_id):
        """Helper function that retrieves a molecule's transitions.

        Args:
            session: SQLAlchemy Session object.
            molecule_id: MoleculeTable integer primary key.

        Returns:
            List of TransitionTable objects.

        Raises:
            TransitionsNotFoundError if no transitions are found.
        """
        stmt = select(TransitionTable).filter_by(molecule_id=molecule_id)
        result = [x[0] for x in session.execute(stmt).all()]
        if not result:
            raise TransitionsNotFoundError(
                      f"transitions not found for molecule {molecule_id}."
                  )
        return result

    def molecules(self):
        """Queries the database for all molecules.

        Returns:
            List of string moelcule chemical formulae.
        """
        with Session(self.engine, future=True) as session:
            stmt = select(MoleculeTable.ordinary_formula)
            return [x[0] for x in session.execute(stmt).all()]

    def gas(self, name):
        """Queries the database for all parameters needed to run a line-by-line calculation.

        Args:
            name: String molecule alias.

        Returns:
            formula: String chemical formula.
            mass: List of float isotopologue masses.
            transitions: List of TransitionTable objects.
            TotalPartionFunction: TotalPartitionFunction object.
        """
        with Session(self.engine, future=True) as session:
            id = self._molecule_id(session, name)
            formula = self._formula(session, id)
            mass = self._mass(session, id)
            transitions = self._transitions(session, id)
            return formula, mass, transitions, TotalPartitionFunction(name, *self.tips(name))

    def tips(self, name):
        """Queries the database for all parameters needed to run TIPS.

        Args:
            name: String molecule alias.

        Returns:
            temperature: Numpy array of temperatures.
            data: Numpy array of data values.

        Raises:
            TipsDataNotFoundError if no TIPS data is found.
        """
        with Session(self.engine, future=True) as session:
            id = self._molecule_id(session, name)
            stmt = select(TipsTable).filter_by(molecule_id=id)
            result = [x[0] for x in session.execute(stmt).all()]
            if not result:
                raise TipsDataNotFoundError(f"no tips data for {name}.")
            data, temperature = [], []
            for value in result:
                data.append(value.data)
                if value.temperature not in temperature:
                    temperature.append(value.temperature)
            data = reshape(asarray(data), (len(data)//len(temperature), len(temperature)))
            temperature = asarray(temperature)
            return temperature, data

    def arts_crossfit(self, name):
        """Queries the database for all parameters needed to run ARTS Crossfit.

        Args:
            name: String molecule alias.

        Returns:
            Path to the cross section dataset.

        Raises:
            CrossSectionNotFoundError if no cross sections are found.
        """
        with Session(self.engine, future=True) as session:
            id = self._molecule_id(session, name)
            stmt = select(ArtsCrossFitTable).filter_by(molecule_id=id)
            try:
                return session.execute(stmt).all()[0][0].path
            except IndexError:
                raise CrossSectionNotFoundError(f"No cross sections for {name}.")


class MoleculeTable(Base):
    """Molecule database table schema."""
    __tablename__ = "molecule"
    id = Column("id", Integer, primary_key=True)
    stoichiometric_formula = Column("stoichiometric_formula", String)
    ordinary_formula = Column("ordinary_formula", String)
    common_name = Column("common_name", String)


class IsotopologueTable(Base):
    """Isotopologue database table schema."""
    __tablename__ = "isotopologue"
    id = Column("id", Integer, primary_key=True)
    molecule_id = Column("molecule_id", Integer, ForeignKey(MoleculeTable.id))
    isoid = Column("isoid", Integer)
    iso_name = Column("iso_name", String)
    abundance = Column("abundance", Float)
    mass = Column("mass", Float)


class MoleculeAliasTable(Base):
    """Molecule alias database table schema."""
    __tablename__ = "molecule_alias"
    id = Column("id", Integer, primary_key=True, autoincrement=True)
    alias = Column("alias", String)
    molecule = Column("molecule", Integer, ForeignKey(MoleculeTable.id))


class TransitionTable(Base):
    """Transition database table schema."""
    __tablename__ = "transition"
    id = Column("id", Integer, primary_key=True, autoincrement=True)
    global_iso_id = Column("global_iso_id", Integer)
    molecule_id = Column("molecule_id", Integer, ForeignKey(MoleculeTable.id))
    local_iso_id = Column("local_iso_id", Integer)
    nu = Column("nu", Float)
    sw = Column("sw", Float)
    gamma_air = Column("gamma_air", Float)
    gamma_self = Column("gamma_self", Float)
    n_air = Column("n_air", Float)
    delta_air = Column("delta_air", Float)
    elower = Column("elower", Float)


class TipsTable(Base):
    """TIPS data table schema."""
    __tablename__ = "tips"
    id = Column("id", Integer, primary_key=True, autoincrement=True)
    molecule_id = Column("molecule_id", Integer, ForeignKey(MoleculeTable.id))
    isotopologue_id = Column("isotopologue_id", Integer)
    temperature = Column("temperature", Float)
    data = Column("data", Float)


class ArtsCrossFitTable(Base):
    """Arts-crossfit table schema."""
    __tablename__ = "artscrossfit"
    id = Column("id", Integer, primary_key=True, autoincrement=True)
    molecule_id = Column("molcule_id", Integer, ForeignKey(MoleculeTable.id))
    path = Column("path", String)


class MetadataTable(Base):
    """Table that describes when data was downloaded."""
    __tablename__ = "metadata"
    id = Column("id", Integer, primary_key=True, autoincrement=True)
    molecule_id = Column("molecule_id", Integer, ForeignKey(MoleculeTable.id))
    database = Column("database", String)
    time = Column("time", String)


class AliasNotFoundError(BaseException):
    pass


class TipsDataNotFoundError(BaseException):
    pass


class IsotopologuesNotFoundError(BaseException):
    pass


class TransitionsNotFoundError(BaseException):
    pass


class CrossSectionNotFoundError(BaseException):
    pass
