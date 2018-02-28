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


def AssignBondOrdersFromTemplateWithoutSanitize(refmol: Mol, mol: Mol) -> Mol:
    """assigns bond orders to a molecule based on the
      bond orders in a template molecule

    AllChem.AssignBondOrdersFromTemplate, performs sanitize which is unwanted,
    so this is a copy with the sanitize removed.
    Also the logging message was reduced from warning to info.

    Args:
        refmol: the template molecule
        mol: the molecule to assign bond orders to

    Returns:
        Molecule with its bond orders overwritten by the template
    """
    refmol2 = Mol(refmol)
    mol2 = Mol(mol)
    # do the molecules match already?
    matching = mol2.GetSubstructMatch(refmol2)
    if not matching:  # no, they don't match
        # check if bonds of mol are SINGLE
        for b in mol2.GetBonds():
            if b.GetBondType() != BondType.SINGLE:
                b.SetBondType(BondType.SINGLE)
                b.SetIsAromatic(False)
        # set the bonds of mol to SINGLE
        for b in refmol2.GetBonds():
            b.SetBondType(BondType.SINGLE)
            b.SetIsAromatic(False)
        # set atom charges to zero;
        for a in refmol2.GetAtoms():
            a.SetFormalCharge(0)
        for a in mol2.GetAtoms():
            a.SetFormalCharge(0)

        matching = mol2.GetSubstructMatches(refmol2, uniquify=False)
        # do the molecules match now?
        if matching:
            if len(matching) > 1:
                logger.info("More than one matching pattern found - picking one")
            matching = matching[0]
            # apply matching: set bond properties
            for b in refmol.GetBonds():
                atom1 = matching[b.GetBeginAtomIdx()]
                atom2 = matching[b.GetEndAtomIdx()]
                b2 = mol2.GetBondBetweenAtoms(atom1, atom2)
                b2.SetBondType(b.GetBondType())
                b2.SetIsAromatic(b.GetIsAromatic())
            # apply matching: set atom properties
            for a in refmol.GetAtoms():
                a2 = mol2.GetAtomWithIdx(matching[a.GetIdx()])
                a2.SetHybridization(a.GetHybridization())
                a2.SetIsAromatic(a.GetIsAromatic())
                a2.SetNumExplicitHs(a.GetNumExplicitHs())
                a2.SetFormalCharge(a.GetFormalCharge())
            if hasattr(mol2, '__sssAtoms'):
                mol2.__sssAtoms = None  # we don't want all bonds highlighted
        else:
            raise ValueError("No matching found")
    return mol2


def hetpdb2mol(pdb_mol: Molecule) -> Mol:
    """Converts atomium Molecule to RDKit molecule via Open Babel

    Via Open Babel because reading PDB with RDKIT looses aromaticity.

    1. Convert PDB with Open Babel to SDF
    2. Read SDF with RDKIT
    3. Assign bond orders based using Open Babel molecule as template

    Args:
        pdb_mol: atomium molecule

    Returns:
        RDKit molecule
    """
    # 1. Convert PDB with Open Babel to SDF
    pdb_block = pdb_mol.to_file_string('pdb')
    babel_mol = pybel.readstring('pdb', pdb_block)
    mol_block = babel_mol.write('mol')
    # 2. Read SDF with RDKIT
    mol_from_babel = MolFromMolBlock(mol_block, removeHs=False)
    # 3. Assign bond orders based using Open Babel molecule as template
    mol_from_rdkit = MolFromPDBBlock(pdb_block, sanitize=False)
    mol_from_rdkit = remove_nonpdb_bonds(mol_from_rdkit, pdb_mol)
    return AssignBondOrdersFromTemplateWithoutSanitize(mol_from_babel, mol_from_rdkit)


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
        reactant = self.rdkit_mol()
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

    def rdkit_mol(self) -> Mol:
        """Return RDKit molecule based on Atomium molecule
        """
        try:
            mol = hetpdb2mol(self.molecule)
        except ValueError as e:
            raise AtomiumParseError(*e.args)
        if not mol:
            raise RdkitParseError('RDKit unable to read ligand ' + self.name())
        try:
            SanitizeMol(mol)
        except ValueError as e:
            # TODO try to fix reactant so it passes SanitizeMol
            raise RdkitParseError(*e.args)
        return mol

    def as_fragment(self):
        return Fragment(self.molecule, self.rdkit_mol())
