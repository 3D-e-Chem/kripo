from math import sqrt
import logging

from atomium.structures import Residue, Atom
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict

from .vector import vector_rotate, normalize, cross_product


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
    axis = (
        pos2[0] - axis_pos[0],
        pos2[1] - axis_pos[1],
        pos2[2] - axis_pos[2],
    )
    d = (
        pos1[0] - pos2[0],
        pos1[1] - pos2[1],
        pos1[2] - pos2[2],
    )
    try:
        dp = vector_rotate(d, angle, axis)
    except:
        dp = d
    dist = sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2])
    return (
        pos2[0] + (dp[0] / dist) * offset,
        pos2[1] + (dp[1] / dist) * offset,
        pos2[2] + (dp[2] / dist) * offset,
    )


def bonded_hydrogens(atom):
    return [a for a in atom.bonded_atoms() if a.element() == 'H']


def bonded_carbons(atom):
    return [a for a in atom.bonded_atoms() if a.element() == 'C']


def sidechain_nitrogens(residue):
    return [a for a in residue.atoms(element='N') if a.name() != 'N']


def sidechain_carbons(residue):
    backbone_carbons = {'C', 'CA'}
    return [a for a in residue.atoms(element='C') if a.name() not in backbone_carbons]


def atoms_by_name(residue, names):
    ring_atoms = set()
    for atom in residue.atoms():
        if atom.name() in names:
            ring_atoms.add(atom)
    return ring_atoms


def center_of_atoms_by_name(residue, names):
    ring_atoms = atoms_by_name(residue, names)
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

    for s_atom in res.atoms(element='S'):
        s = s_atom.location()
        (s1_atom, s2_atom) = list(s_atom.bonded_atoms())
        s1 = s1_atom.location()
        s2 = s2_atom.location()
        s1_s = (
            s1[0] - s[0],
            s1[1] - s[1],
            s1[2] - s[2],
        )
        s2_s = (
            s2[0] - s[0],
            s2[1] - s[1],
            s2[2] - s[2],
        )
        logging.warning(
            "Adding hydrogens to sulfur with 90 degree angle between h-s-h, instead of 109.471 degrees")
        anti_s12 = (
            (s1_s[0] + s2_s[0]) * -0.5,
            (s1_s[1] + s2_s[1]) * -0.5,
            (s1_s[2] + s2_s[2]) * -0.5,
        )
        plane = normalize(
            cross_product(s1_s, s2_s)
        )
        plane_dist = (
            sqrt(3 / 4) * plane[0],
            sqrt(3 / 4) * plane[1],
            sqrt(3 / 4) * plane[2],
        )
        h2 = (
            s[0] - plane_dist[0] + anti_s12[0],
            s[1] - plane_dist[1] + anti_s12[1],
            s[2] - plane_dist[2] + anti_s12[2],
        )
        h3 = (
            s[0] + plane_dist[0] + anti_s12[0],
            s[1] + plane_dist[1] + anti_s12[1],
            s[2] + plane_dist[2] + anti_s12[2],
        )
        atom_h2 = Atom(element='H',
                       x=h2[0],
                       y=h2[1],
                       z=h2[2],
                       name='H')
        res.add_atom(atom_h2)
        s_atom.bond(atom_h2)
        atom_h3 = Atom(element='H',
                       x=h3[0],
                       y=h3[1],
                       z=h3[2],
                       name='H')
        res.add_atom(atom_h3)
        s_atom.bond(atom_h3)
    return res
