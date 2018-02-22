import logging
import pkg_resources
from typing import Dict

from rdkit.Chem import AllChem, MolToSmiles, SanitizeMol, Mol, Atom, MolToMolBlock
from rdkit.Geometry.rdGeometry import Point3D


def map_atoms(mol_a: Mol, mol_b: Mol) -> Dict[int, int]:
    """Find all atoms in mol_a in mol_b

    Args:
        mol_a:
        mol_b:

    Returns:
        Dictionary where key are atom indices from mol_a and values are atom indices from mol_b
    """
    a_conf = mol_a.GetConformer()
    b_conf = mol_b.GetConformer()

    # find mapping of atoms from mol to parent
    a2b_atoms = {}
    b2a_atoms = {}
    for a in mol_a.GetAtoms():
        a_pos = a_conf.GetAtomPosition(a.GetIdx())
        for b in mol_b.GetAtoms():
            if b.GetIdx() in b2a_atoms:
                # Don't try to map same atom twice
                continue
            bp_pos = b_conf.GetAtomPosition(b.GetIdx())
            if b.Match(a) and a_pos.x == bp_pos.x and a_pos.y == bp_pos.y and a_pos.z == bp_pos.z:
                a2b_atoms[a.GetIdx()] = b.GetIdx()
                b2a_atoms[b.GetIdx()] = a.GetIdx()

    return a2b_atoms


def embed_r_groups(mol: Mol, parent: Mol):
    """Assign coordinates to R groups

    Args:
        mol: Molecule with atoms with '*' as symbol
        parent: Molecule that was used to generate `mol`

    Raises:
        LookupError: When atoms bonded to R group can not be found in parent

    """
    mol2parent_idxs = map_atoms(mol, parent)
    parent2mol_idxs = {v for v in mol2parent_idxs.values()}
    idx2parent_atom = {a.GetIdx(): a for a in parent.GetAtoms()}

    mol_conf = mol.GetConformer()
    parent_conf = parent.GetConformer()
    r_pos_cache = set()
    for r_atom in mol.GetAtoms():
        if r_atom.GetSymbol() == '*':
            for r_bond in r_atom.GetBonds():
                # atom bonded to r group in mol
                r_bonded_atom = r_bond.GetOtherAtom(r_atom)
                r_bonded_idx = r_bonded_atom.GetIdx()
                if r_bonded_idx in mol2parent_idxs:
                    # same atom in parent
                    r_bonded_atom_parent_idx = mol2parent_idxs[r_bonded_idx]
                    r_bonded_atom_parent = idx2parent_atom[r_bonded_atom_parent_idx]
                    for r_bonded_atom_parent_bond in r_bonded_atom_parent.GetBonds():
                        # Atom bonded to atom which is bonded to R group in parent
                        r_bonded_atom_parend_bonded_atom = r_bonded_atom_parent_bond.GetOtherAtom(r_bonded_atom_parent)
                        r_bonded_atom_parend_bonded_atom_idx = r_bonded_atom_parend_bonded_atom.GetIdx()
                        if r_bonded_atom_parend_bonded_atom_idx not in parent2mol_idxs:
                            # atom in parent but not in mol which was replaced R group
                            point = parent_conf.GetAtomPosition(r_bonded_atom_parend_bonded_atom_idx)
                            if serialize_point(point) not in r_pos_cache:
                                # Position has not been used for other R group
                                mol_conf.SetAtomPosition(r_atom.GetIdx(), point)
                                r_pos_cache.add(serialize_point(point))
                                break
                else:
                    raise LookupError('Atom bonded to R group not found in parent')


def serialize_point(point: Point3D):
    return str(point.x) + str(point.y) + str(point.z)


class Reactor:
    """Reactor, facilitator of reactions

    Attributes:
        steps: Maximum number of times a product is reacted

    """
    def __init__(self, steps=99):
        self.reactions = []
        self.steps = steps
        self.load_reactions()

    def load_reactions(self):
        stream = pkg_resources.resource_stream('kripo', 'data/combined.smirks')
        for line in stream:
            cols = line.decode('ascii').strip().split()
            if len(cols) > 0 and cols[0]:
                smirk = cols[0]
                reaction = AllChem.ReactionFromSmarts(smirk)
                self.reactions.append(reaction)

    def react(self, reactant):
        products = set()
        product_smiles = set()
        n = self.steps
        new_mols = [reactant]
        while n > 0 and new_mols != []:
            mols = new_mols
            new_mols = []
            for mol in mols:
                SanitizeMol(mol)
                for reaction in self.reactions:
                    for ps in reaction.RunReactants((mol,)):
                        q = ps[0]
                        SanitizeMol(q)
                        smile = MolToSmiles(q)
                        if smile not in product_smiles:
                            embed_r_groups(q, mol)
                            new_mols.append(q)
                            product_smiles.add(smile)
                            products.add(q)

            n -= 1

        return products
