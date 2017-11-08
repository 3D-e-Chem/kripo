from hashlib import md5

from atomium.structures.chains import Site
from atomium.structures.molecules import Molecule
from rdkit.Chem import Mol, MolToSmiles, MolToMolBlock

"""Residues within radius of ligand are site residues"""
BINDING_SITE_RADIUS = 6


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

    def atoms(self):
        """Atoms of fragment

        Yields:
            atoms.html#atomium.structures.atoms.Atom: An atom
        """
        fragment_names = set(self.atom_names())
        for atom in self.parent.atoms(exclude="H"):
            if atom.name() in fragment_names:
                yield atom

    def site(self, radius=BINDING_SITE_RADIUS) -> Site:
        """Site of fragment

        Args:
            radius (float): Radius of ligand within residues are included in site

        Returns:
            atomium.structures.chains.Site: Site
        """
        atoms = self.parent.atoms(exclude="H")
        nearby = set()
        # assume atom names of parent have been retained in fragment rdkit molecule
        fragment_names = set(self.atom_names())
        for atom in atoms:
            if atom.name() in fragment_names:
                nearby.update(atom.nearby(radius, exclude="H"))
        # residues excluding parent atoms
        residues = [atom.residue() for atom in nearby if atom not in atoms]

        # TODO add hydrogens?
        return Site(*residues, ligand=self.parent)

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
