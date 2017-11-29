from math import sqrt
import logging

from atomium.structures import Residue
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict

from ..protonate import protonate_protein, protonate_ligand
from .vector import vector_rotate


def feature_pos_of_bond(atom1, atom2, offset):
    pos1 = atom1.location()
    pos2 = atom2.location()
    d = (
        pos1[0] - pos2[0],
        pos1[1] - pos2[1],
        pos1[2] - pos2[2],
    )
    dist = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
    return (
        pos1[0] + (d[0] / dist) * offset,
        pos1[1] + (d[1] / dist) * offset,
        pos1[2] + (d[2] / dist) * offset,
    )


def feature_pos_of_bond_rotated(atom1, atom2, offset, angle, axis_atom):
    pos1 = atom1.location()
    pos2 = atom2.location()
    axis_pos = axis_atom.location()
    d = (
        axis_pos[0] - pos2[0],
        axis_pos[1] - pos2[1],
        axis_pos[2] - pos2[2],
    )
    try:
        dp = vector_rotate(d, angle, d)
    except:
        dp = d
    dist = sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2])
    return (
        pos1[0] + (dp[0] / dist) * offset,
        pos1[1] + (dp[1] / dist) * offset,
        pos1[2] + (dp[2] / dist) * offset,
    )


def bonded_hydrogens(atom):
    return [a for a in atom.bonded_atoms() if a.element() == 'H']


def sidechain_nitrogens(residue):
    return [a for a in residue.atoms(element='N') if a.name() != 'N']


def sidechain_carbons(residue):
    backbone_carbons = {'C', 'CA'}
    return [a for a in residue.atoms(element='C') if a.name() not in backbone_carbons]


def atoms_of_ring(residue, names):
    ring_atoms = set()
    for atom in residue.atoms():
        if atom.name() in names:
            ring_atoms.add(atom)
    return ring_atoms


def center_of_ring(residue, names):
    ring_atoms = atoms_of_ring(residue, names)
    ring_center = (
        sum([a.location()[0] for a in ring_atoms]) / len(ring_atoms),
        sum([a.location()[1] for a in ring_atoms]) / len(ring_atoms),
        sum([a.location()[2] for a in ring_atoms]) / len(ring_atoms),
    )
    return ring_center


def copy_residue(residue: Residue) -> Residue:
    block = residue.to_file_string('pdb')
    pdb = pdb_dict_to_pdb(pdb_string_to_pdb_dict(block))
    return pdb.model().residue()


def add_hydrogens2sulfur_as_carbon(resin: Residue) -> Residue:
    res = copy_residue(resin)

    for a in res.atoms(element='S'):
        a.element('C')
        a.name(a.name().replace('S', 'C'))
        # for h in bonded_hydrogens(a):
        #     res.remove_atom(h)
    logging.warning("TODO: Not adding hydrogens to sulfur, less features will be generated until this is implemented")
    return res

    # TODO find way to add hydrogens to carbon that was prev a sulfur

    max_serial_number = max([a.atom_id() for a in res.atoms() if a.atom_id()])
    h_atom_id = max_serial_number + 1
    for a in res.atoms():
        if not a.atom_id():
            a._id = h_atom_id
            h_atom_id += 1

    block = res.to_file_string('pdb')
    """
    * reduce: does not add Hs
    * openbabel: ignores existing hs and re-adds all of them
    * rdkit: ignores bonds when loading pdb block
    
    """
    hblock = protonate_ligand(block)
    hpdb = pdb_dict_to_pdb(pdb_string_to_pdb_dict(hblock))
    return hpdb.model().residue()
