from sqlalchemy import Column, create_engine, Float, ForeignKey, Integer, select, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import Session

from .webapi import NoIsotopologueError, NoTransitionsError


Base = declarative_base()


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
        self.engine = create_engine(f"sqlite+pysqlite:///{path}", echo=echo, future=True)
        Base.metadata.create_all(self.engine)

    def create(self, webapi, molecules="all"):
        """Downloads data from HITRAN and inserts it into the database tables.

        Args:
            webapi: WebApi object.
            molecules: List of string molecule chemical formulae.
        """
        with Session(self.engine, future=True) as session:
            all_molecules = webapi.download_molecules()
            for i, molecule in enumerate(all_molecules):

                # Write out progress.
                print("Working on molecule {} / {} ({})".format(
                    i + 1, len(all_molecules), molecule.ordinary_formula)
                )

                # Support for only using a subset of molecules.
                if molecules != "all":
                    if molecule.ordinary_formula not in molecules:
                        continue

                # Store the molecule metadata.
                session.add(
                    MoleculeTable(
                        id=molecule.id,
                        stoichiometric_formula=molecule.stoichiometric_formula,
                        ordinary_formula=molecule.ordinary_formula,
                        common_name=molecule.common_name
                    )
                )

                # Store the molecule aliases.
                aliases = [x["alias"] for x in molecule.aliases]
                for alias in aliases:
                    session.add(
                        MoleculeAliasTable(
                            alias=alias,
                            molecule=molecule.id
                        )
                    )

                # Store isotopologues.
                isotopologues = webapi.download_isotopologues(molecule)
                for isotopologue in isotopologues:
                    session.add(
                        IsotopologueTable(
                            id=isotopologue.id,
                            molecule_id=molecule.id,
                            isoid=isotopologue.isoid,
                            iso_name=isotopologue.iso_name,
                            abundance=isotopologue.abundance,
                            mass=isotopologue.mass
                        )
                    )

                # Store the transitions.
                parameters = ["global_iso_id", "molec_id", "local_iso_id", "nu", "sw",
                              "gamma_air", "gamma_self", "n_air", "delta_air", "elower"]
                try:
                    transitions = webapi.download_transitions(isotopologues, 0., 1.e8, parameters)
                except NoIsotopologueError:
                    print("No isotopologues for molecule {}.".format(molecule.ordinary_formula))
                    continue
                except NoTransitionsError:
                    print("No transitions for molecule {}.".format(molecule.ordinary_formula))
                    continue
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
                            elower=transition.elower
                        )
                    )
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

    def _mass(self, session, molecule_id):
        """Helper function that retrieves a molecule's isotopologue masses.

        Args:
            session: SQLAlchemy Session object.
            molecule_id: MoleculeTable integer primary key.

        Returns:
            List of float isotopologue masses.
        """
        stmt = select(IsotopologueTable.mass).filter_by(molecule_id=molecule_id)
        return [x[0] for x in session.execute(stmt).all()]

    def _molecule_id(self, session, name):
        """Helper function that retrieves a molecules id.

        Args:
            session: SQLAlchemy Session object.
            name: String molecule alias.

        Returns:
            MoleculeTable integer primary key.
        """
        stmt = select(MoleculeAliasTable.molecule).filter_by(alias=name)
        return session.execute(stmt).all()[0][0]

    def _transitions(self, session, molecule_id):
        """Helper function that retrieves a molecule's transitions.

        Args:
            session: SQLAlchemy Session object.
            molecule_id: MoleculeTable integer primary key.

        Returns:
            List of TransitionTable objects.
        """
        stmt = select(TransitionTable).filter_by(molecule_id=molecule_id)
        return [x[0] for x in session.execute(stmt).all()]

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
        """
        with Session(self.engine, future=True) as session:
            id = self._molecule_id(session, name)
            formula = self._formula(session, id)
            mass = self._mass(session, id)
            transitions = self._transitions(session, id)
            return formula, mass, transitions


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
