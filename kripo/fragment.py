from hashlib import md5
from typing import Set, List, Tuple

from atomium.structures.chains import Site
from atomium.structures.molecules import Molecule, Residue
from atomium.structures.atoms import Atom
from math import sqrt
from rdkit.Chem import Mol, MolToSmiles, RemoveHs

"""Residues within radius of ligand are site residues"""
BINDING_SITE_RADIUS = 6


def distance_between_positions(pos1, pos2):
    d = (
        pos1[0] - pos2[0],
        pos1[1] - pos2[1],
        pos1[2] - pos2[2],
    )
    return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])


def is_residue_nearby(fragment_positions: List[Tuple[float, float, float]], residue: Residue, radius: float) -> bool:
    residue_atoms = residue.atoms()
    min_distance = 9999.0
    if radius > min_distance:
        raise ValueError("Radius must be smaller than {0}".format(min_distance))
    for fragment_position in fragment_positions:
        for residue_atom in residue_atoms:
                dist = distance_between_positions(fragment_position, residue_atom.location())
                if dist < min_distance:
                    min_distance = dist
    return min_distance < radius


class Fragment:
    """Fragment of a ligand

    Attributes:
        parent (atomium.structures.molecules.Molecule): The parent ligand
        molecule (Mol): Fragment molecule with hydrogens

    """
    def __init__(self, parent: Molecule, molecule: Mol):
        self.parent = parent
        self.molecule = molecule

    def atom_names(self, include_hydrogen=True):
        """Ligand atom names which make up the fragment

        Excludes hydrogens

        Returns:
            List[str]: Atom names

        """
        raise NotImplemented('Mapping of fragment atoms to pdb ligand atoms')

    def atoms(self) -> Set[Atom]:
        """Atoms of fragment

        Returns:
            collection of atoms
        """
        raise NotImplemented('Mapping of fragment atoms to pdb ligand atoms')

    def site(self, radius=BINDING_SITE_RADIUS) -> Site:
        """Site of fragment

        If part of residue is inside radius then it is included in site.

        Args:
            radius (float): Radius of ligand within residues are included in site

        Returns:
            atomium.structures.chains.Site: Site
        """
        fragment_positions = self.molecule.GetConformer().GetPositions()
        atoms_of_near_residues = set()
        residues = self.parent.model().residues()
        for residue in residues:
            if is_residue_nearby(fragment_positions, residue, radius):
                atoms_of_near_residues.update(residue.atoms())

        return Site(*atoms_of_near_residues)

    def nr_r_groups(self):
        """Number of R groups in fragment

        Returns:
            int: number of R groups

        """
        counter = 0
        for atom in self.molecule.GetAtoms():
            if atom.GetSymbol() == '*':
                counter += 1

        return counter

    def smiles(self):
        """Smiles string of fragment

        Returns:
            str: Smiles
        """
        return MolToSmiles(self.molecule)

    def hash_code(self):
        """Hash code of fragment

        Returns:
            str: hash code
        """
        smiles = self.smiles().encode('ascii')
        return md5(smiles).hexdigest()

    @property
    def name(self) -> str:
        """Name of fragment"""
        try:
            return self.molecule.GetProp('_Name')
        except KeyError:
            return ""

    @name.setter
    def name(self, name: str):
        self.molecule.SetProp('_Name', name)

    def unprotonated_molecule(self) -> Mol:
        """Return molecule with all hydrogens removed
        """
        return RemoveHs(self.molecule)
