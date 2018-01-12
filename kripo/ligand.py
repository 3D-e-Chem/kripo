from typing import List

from atomium.structures.chains import Site
from atomium.structures.molecules import Molecule
from rdkit.Chem import MolFromPDBBlock, SanitizeMol
from rdkit.Chem.Descriptors import MolWt

from .reactor import Reactor
from .fragment import Fragment, BINDING_SITE_RADIUS


class RdkitParseError(ValueError):
    pass


class Ligand:
    """Ligand of a ligand-protein complex

    Attributes:
        molecule (atomium.structures.molecules.Molecule): Atomium molecule

    """
    def __init__(self, molecule: Molecule):
        self.molecule = molecule

    def name(self):
        """Hetero code of the ligand"""
        return self.molecule.name()

    def site(self, radius=BINDING_SITE_RADIUS):
        """Site of ligand

        Args:
            radius (float): Radius of ligand within residues are included in site
        Returns:
            atomium.structures.chains.Site: Site

        """
        atoms, nearby = self.molecule.atoms(exclude="H"), set()
        for atom in atoms:
            nearby.update(atom.nearby(radius, exclude="H"))
        residues = [atom.residue() for atom in nearby if atom not in atoms]
        residues = [residue for residue in residues if residue]
        return Site(*residues, ligand=self.molecule)

    def pdb_block(self):
        """Ligand as PDB block

        Returns:
            str: PDB block
        """
        return self.molecule.to_file_string('pdb')

    def fragments(self) -> List[Fragment]:
        """Fragments ligand using RDKit.

        Includes self as first element.

        Returns:
            List[Fragment]: Ordered by weight, heaviest first

        """
        block = self.pdb_block()
        reactant = MolFromPDBBlock(block, sanitize=False)
        if not reactant:
            raise RdkitParseError('RDKit unable to read ligand ' + self.name())
        try:
            SanitizeMol(reactant)
        except ValueError as e:
            # TODO try to fix reactant so it passes SanitizeMol
            raise RdkitParseError(*e.args)
        mols = [reactant]
        products = Reactor().react(reactant)
        mols.extend(products)
        mols.sort(key=lambda m: MolWt(m), reverse=True)
        return [Fragment(self.molecule, mol) for mol in mols]

    def id(self):
        """Identifier of ligand

        Returns:
            str: Unique id of ligand
        """
        return self.molecule.molecule_id()

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
