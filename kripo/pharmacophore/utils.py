from math import sqrt

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
    return [hyd for hyd in atom.bonded_atoms() if hyd.element() == 'H']


def atoms_of_ring(residue, names):
    ring_atoms = set()
    for atom in residue.atoms():  # TODO only sidechain
        if atom.name() in names:
            continue
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
