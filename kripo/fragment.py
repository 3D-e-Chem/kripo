from hashlib import md5
from typing import List, Set, Tuple

from atomium.structures.chains import Site
from atomium.structures.molecules import Molecule, Residue
from atomium.structures.atoms import Atom
from math import sqrt
from rdkit.Chem import Mol, MolToSmiles, RemoveHs

"""Residues within radius of ligand are site residues"""
BINDING_SITE_RADIUS = 6


def bounding_box(atoms: List[Tuple[float, float, float]]):
    min_x = float('inf')
    max_x = float('-inf')
    min_y = float('inf')
    max_y = float('-inf')
    min_z = float('inf')
    max_z = float('-inf')
    if not atoms:
        raise ValueError('List can not be empty')
    for atom in atoms:
        if atom[0] < min_x:
            min_x = atom[0]
        if atom[1] < min_y:
            min_y = atom[1]
        if atom[2] < min_z:
            min_z = atom[2]
        if atom[0] > max_x:
            max_x = atom[0]
        if atom[1] > max_y:
            max_y = atom[1]
        if atom[2] > max_z:
            max_z = atom[2]

    return (min_x, max_x), (min_y, max_y), (min_z, max_z)


def bounding_boxes_overlap(fragment_bounding_box, radius, residue_bounding_box):
    return fragment_bounding_box[0][0] - radius <= residue_bounding_box[0][1] and \
           fragment_bounding_box[0][1] + radius >= residue_bounding_box[0][0] and \
           fragment_bounding_box[1][0] - radius <= residue_bounding_box[1][1] and \
           fragment_bounding_box[1][1] + radius >= residue_bounding_box[1][0] and \
           fragment_bounding_box[2][0] - radius <= residue_bounding_box[2][1] and \
           fragment_bounding_box[2][1] + radius >= residue_bounding_box[2][0]


def distance_between_positions(pos1, pos2):
    d = (
        pos1[0] - pos2[0],
        pos1[1] - pos2[1],
        pos1[2] - pos2[2],
    )
    return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])


def is_residue_nearby(fragment_positions: List[Tuple[float, float, float]], residue: Residue, radius: float) -> bool:
    residue_positions = [a.location() for a in residue.atoms()]
    residue_bounding_box = bounding_box(residue_positions)
    fragment_bounding_box = bounding_box(fragment_positions)
    if not bounding_boxes_overlap(fragment_bounding_box, radius, residue_bounding_box):
        return False
    square_radius = radius * radius
    for fragment_atom_position in fragment_positions:
        for residue_atom_position in residue_positions:
            x = residue_atom_position[0] - fragment_atom_position[0]
            y = residue_atom_position[1] - fragment_atom_position[1]
            z = residue_atom_position[2] - fragment_atom_position[2]
            dist = x * x + y * y + z * z
            if dist < square_radius:
                return True
    return False


class Fragment:
    """Fragment of a ligand

    Attributes:
        parent (atomium.structures.molecules.Molecule): The parent ligand
        molecule (Mol): Fragment molecule with hydrogens

    """
    def __init__(self, parent: Molecule, molecule: Mol):
        self.parent = parent
        self.molecule = molecule

    def atom_names(self):
        """Ligand atom names which make up the fragment

        Excludes hydrogens and R groups

        Returns:
            List[str]: Atom names

        """
        return [a.name() for a in self.atoms()]

    def atoms(self) -> Set[Atom]:
        """Atoms of fragment

        Excludes hydrogens and R groups

        Returns:
            collection of atoms
        """
        theset = set()
        conf = self.molecule.GetConformer()
        ignored_symbols = {'H', '*'}
        delta = 0.002
        for a in self.molecule.GetAtoms():
            if a.GetSymbol() in ignored_symbols:
                continue

            frag_symbol = a.GetSymbol().upper()
            p = conf.GetAtomPosition(a.GetIdx())
            frag_loc = (p.x, p.y, p.z)

            # Find fragment atom in whole ligand aka parent based on symbol and position
            for parent_atom in self.parent.atoms():
                if parent_atom in theset:
                    continue
                if parent_atom.element() != frag_symbol:
                    continue
                parent_loc = parent_atom.location()
                if all((
                    abs(parent_loc[0] - frag_loc[0]) < delta,
                    abs(parent_loc[1] - frag_loc[1]) < delta,
                    abs(parent_loc[2] - frag_loc[2]) < delta,
                )):
                    # Match found
                    theset.add(parent_atom)
                    break

        return theset

    def atom_positions(self) -> List[Tuple[float, float, float]]:
        conf = self.molecule.GetConformer()
        # Using conf.GetPositions gave segfault, so iterate over each atom
        return [(p.x, p.y, p.z) for p in [conf.GetAtomPosition(i) for i in range(conf.GetNumAtoms())]]

    def site(self, radius=BINDING_SITE_RADIUS) -> Site:
        """Site of fragment

        If part of residue is inside radius then it is included in site.

        Args:
            radius (float): Radius of ligand within residues are included in site

        Returns:
            atomium.structures.chains.Site: Site
        """
        fragment_positions = self.atom_positions()
        atoms_of_near_residues = set()
        residues = self.parent.model().residues()
        for residue in residues:
            if is_residue_nearby(fragment_positions, residue, radius):
                atoms_of_near_residues.update(residue.atoms())

        name = 'fragment of ' + self.parent.name()
        ligand = Molecule(*self.atoms(), name=name)
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
