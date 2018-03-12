import logging
from typing import List

import pybel
from atomium.structures.chains import Site
from atomium.structures.molecules import Molecule
from rdkit.Chem import MolFromPDBBlock, SanitizeMol, RWMol, Mol, MolFromMol2Block, MolFromMolBlock, BondType
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from rdkit.Chem.Descriptors import MolWt

from .reactor import Reactor
from .fragment import Fragment, BINDING_SITE_RADIUS

logger = logging.getLogger(__name__)


class RdkitParseError(ValueError):
    pass


class AtomiumParseError(ValueError):
    pass


class Ligand:
    """Ligand of a ligand-protein complex

    Attributes:
        atomium_mol: Atomium molecule
        rdkit_mol: protonated RDKit molecule

    """
    def __init__(self, atomium_mol: Molecule, rdkit_mol: Mol):
        self.atomium_mol = atomium_mol
        self.rdkit_mol = rdkit_mol

    def name(self):
        """Hetero code of the ligand"""
        return self.atomium_mol.name()

    def fragments(self) -> List[Fragment]:
        """Fragments ligand using RDKit.

        Includes self as first element.

        Returns:
            List[Fragment]: Ordered by weight, heaviest first

        """
        reactant = self.rdkit_mol
        mols = [reactant]
        products = Reactor().react(reactant)
        mols.extend(products)
        mols.sort(key=lambda m: MolWt(m), reverse=True)
        return [Fragment(self.atomium_mol, mol) for mol in mols]

    def id(self):
        """Identifier of ligand

        Returns:
            str: Unique id of ligand
        """
        return self.atomium_mol.molecule_id()

    def chain(self):
        """Chain of ligand

        Returns:
            str: Chain letter

        """
        return self.id()[0]

    def seq_nr(self):
        """Sequence number of ligand

        Returns:
            int: sequence number
        """
        return int(self.id()[1:])

    def as_fragment(self):
        return Fragment(self.atomium_mol, self.rdkit_mol)
