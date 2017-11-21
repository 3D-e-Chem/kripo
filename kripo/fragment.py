from hashlib import md5
from typing import Set

from atomium.structures.chains import Site
from atomium.structures.molecules import Molecule, Residue
from atomium.structures.atoms import Atom
from rdkit.Chem import Mol, MolToSmiles, MolToMolBlock

"""Residues within radius of ligand are site residues"""
BINDING_SITE_RADIUS = 6


def is_residue_nearby(fragment_atoms: Set[Atom], residue: Residue, radius: float) -> bool:
    residue_atoms = residue.atoms()
    min_distance = 9999.0
    if radius > min_distance:
        raise ValueError("Radius must be smaller than {0}".format(min_distance))
    for fragment_atom in fragment_atoms:
        for residue_atom in residue_atoms:
                dist = fragment_atom.distance_to(residue_atom)
                if dist < min_distance:
                    min_distance = dist
    return min_distance < radius


class Fragment:
    """Fragment of a ligand

    Attributes:
        parent (atomium.structures.molecules.Molecule): The parent ligand
        molecule (Mol): Fragment molecule

    """
    def __init__(self, parent: Molecule, molecule: Mol):
        self.parent = parent
        self.molecule = molecule

    def atom_names(self):
        """Ligand atom names which make up the fragment

        Returns:
            List[str]: Atom names

        """
        return [a.GetPDBResidueInfo().GetName().strip() for a in self.molecule.GetAtoms() if a.GetPDBResidueInfo() is not None]

    def atoms(self) -> Set[Atom]:
        """Atoms of fragment

        Returns:
            collection of atoms
        """
        fragment_names = set(self.atom_names())
        atoms = set()

        for atom in self.parent.atoms():
            if atom.name() in fragment_names:
                atoms.add(atom)
        return atoms

    def site(self, radius=BINDING_SITE_RADIUS) -> Site:
        """Site of fragment

        If part of residue is inside radius then it is included in site.

        Args:
            radius (float): Radius of ligand within residues are included in site

        Returns:
            atomium.structures.chains.Site: Site
        """
        fragment_atoms = self.atoms()
        atoms_of_near_residues = set()
        residues = self.parent.model().residues()
        for residue in residues:
            if is_residue_nearby(fragment_atoms, residue, radius):
                atoms_of_near_residues.update(residue.atoms())

        ligand = Molecule(*fragment_atoms, molecule_id=self.parent.molecule_id(), name=self.parent.name())
        return Site(*atoms_of_near_residues, ligand=ligand)

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

    def mol_block(self, name):
        """Mol block of fragment

        Args:
            name (str): Name inserted into mol block

        Returns:
            str: Mol block
        """
        mol = Mol(self.molecule)
        mol.SetProp('_Name', name)
        return MolToMolBlock(mol)
