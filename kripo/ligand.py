import logging
from typing import List

from atomium.structures.chains import Site
from atomium.structures.molecules import Molecule
from rdkit.Chem import MolFromPDBBlock, SanitizeMol, RWMol, Mol
from rdkit.Chem.Descriptors import MolWt

from .reactor import Reactor
from .fragment import Fragment, BINDING_SITE_RADIUS

logger = logging.getLogger()


class RdkitParseError(ValueError):
    pass


class AtomiumParseError(ValueError):
    pass


def remove_nonpdb_bonds(rdkit_mol: Mol, atomium_mol: Molecule) -> Mol:
    """Checks if bonds in RDKit molecule are also present in PDB. If absent then the bond is removed.

    Args:
        rdkit_mol: RDKit molecule to prune
        atomium_mol: Reference PDB molecule

    Returns:
        RDKit molecule which has been pruned
    """
    rwmol = RWMol(rdkit_mol)
    for a in rwmol.GetAtoms():
        a_idx = a.GetIdx()
        a_serial = a.GetPDBResidueInfo().GetSerialNumber()
        a_atom = atomium_mol.atom(a_serial)
        if a.GetExplicitValence() != len(a_atom.bonds()):
            for bond in a.GetBonds():
                o = bond.GetOtherAtom(a)
                o_idx = o.GetIdx()
                o_serial = o.GetPDBResidueInfo().GetSerialNumber()
                o_atom = atomium_mol.atom(o_serial)
                if not a_atom.bond_with(o_atom):
                    rwmol.RemoveBond(a_idx, o_idx)
                    logger.info("Removing bond %s:%s - %s:%s from RDKit molecule as it is not present in PDB",
                                a_serial, a.GetSymbol(),
                                o_serial, o.GetSymbol())
    return rwmol.GetMol()


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
        try:
            block = self.pdb_block()
        except ValueError as e:
            raise AtomiumParseError(*e.args)
        reactant = MolFromPDBBlock(block, sanitize=False)
        if not reactant:
            raise RdkitParseError('RDKit unable to read ligand ' + self.name())
        reactant = remove_nonpdb_bonds(reactant, self.molecule)
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
