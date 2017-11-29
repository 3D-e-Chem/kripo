import logging

from atomium.structures import Residue

from .feature import Feature
from .utils import feature_pos_of_bond, feature_pos_of_bond_rotated, bonded_hydrogens, atoms_of_ring, \
    center_of_ring, sidechain_nitrogens, sidechain_carbons, add_hydrogens2sulfur_as_carbon
from .vector import center_of_triangle, above, below

H_dist = 0.8  # HYDROPHOBE
Rp_width = 3.4  # PLACEMENT OF PI AROMATIC FEATURES
Rp_dist = 4.0  # AROMATIC PI STACKING
Rt_dist = -1  # ATOMATIC T STACKING
A_dist = 0.8  # ACCEPTORS
O_dist = 0.2  # DONORS
P_dist = 0.3  # POSITIVE CHARGES
P_width = 0.0  # POSITIVE CHARGES
N_dist = 1.4  # NEGATIVE CHARGES
N_width = 1.0  # NEGATIVE CHARGES


def distance_to(atom, ligand):
    min_dist = None
    for a in ligand.atoms():
        dist = atom.distance_to(a)
        if min_dist is None or dist < min_dist:
            min_dist = dist
    return min_dist


def closest_to(atoms, ligand):
    closest = None
    min_dist = None
    if len(atoms) == 1:
        return atoms[0]
    for a in atoms:
        dist = distance_to(a, ligand)
        if min_dist is None or dist < min_dist:
            min_dist = dist
            closest = a
    return closest


def features_from_backbone_amine(residue: Residue, ligand):
    nitrogen = residue.atom(name='N')
    hs = bonded_hydrogens(nitrogen)
    hydrogen = closest_to(hs, ligand)
    if not hydrogen:
        return set()

    feature_pos = feature_pos_of_bond(hydrogen, nitrogen, O_dist)
    feature = Feature('HACC', feature_pos)

    return {feature}


def features_from_backbone_carbonyl(residue):
    oxygen = residue.atom(name='O')
    carbon = residue.atom(name='C')

    feature_pos = feature_pos_of_bond(oxygen, carbon, A_dist)
    feature = Feature('HDON', feature_pos)

    return {feature}


def features_from_alanine_sidechain(residue):
    cb = residue.atom(name='CB')
    ca = residue.atom(name='CA')

    feature_pos = feature_pos_of_bond(cb, ca, H_dist)
    feature = Feature('LIPO', feature_pos)
    features = {feature}
    for hyd in bonded_hydrogens(cb):
        feature_pos = feature_pos_of_bond(hyd, cb, H_dist)
        feature = Feature('LIPO', feature_pos)
        features.add(feature)

    return features


def features_from_sidechain_hydrophobes(residue):
    cb = residue.atom(name='CB')
    cg = residue.atom(name='CG')
    features = set()
    for hyd in bonded_hydrogens(cb):
        feature_pos = feature_pos_of_bond(hyd, cb, H_dist)
        feature = Feature('LIPO', feature_pos)
        features.add(feature)
    for hyd in bonded_hydrogens(cg):
        feature_pos = feature_pos_of_bond(hyd, cg, H_dist)
        feature = Feature('LIPO', feature_pos)
        features.add(feature)
    return features


def features_from_sidechain_donors(residue):
    features = set()
    for nitrogen in sidechain_nitrogens(residue):
        for hydrogen in bonded_hydrogens(nitrogen):
            feature_pos = feature_pos_of_bond(hydrogen, nitrogen, O_dist)
            feature = Feature('HACC', feature_pos)
            features.add(feature)
    return features


def features_from_sidechain_positives(residue):
    cz = residue.atom(name='CZ')
    nh1 = residue.atom(name='NH1')
    nh2 = residue.atom(name='NH2')

    cp = center_of_triangle(nh1, cz, nh2)

    features = set()
    for nitrogen in sidechain_nitrogens(residue):
        middle_pos = feature_pos_of_bond(nitrogen, cz, P_dist)
        features |= {
            Feature('POSC', above(middle_pos, cp, P_width)),
            Feature('POSC', below(middle_pos, cp, P_width))
        }

    return features


def features_from_asparagine_sidechain(residue):
    cg = residue.atom(name='CG')
    od1 = residue.atom(name='OD1')
    nd2 = residue.atom(name='ND2')

    feature_pos = feature_pos_of_bond(od1, cg, A_dist)
    feature = Feature('HDON', feature_pos)
    features = {feature}

    for hyd in bonded_hydrogens(nd2):
        feature_pos = feature_pos_of_bond(hyd, nd2, O_dist)
        feature = Feature('HACC', feature_pos)
        features.add(feature)

    return features


def features_from_asparticacids_sidechain(residue):
    cb = residue.atom(name='CB')
    cg = residue.atom(name='CG')
    od1 = residue.atom(name='OD1')
    od2 = residue.atom(name='OD2')

    cp = center_of_triangle(od1, cg, od2)

    features = {Feature('HDON', feature_pos_of_bond(od1, cg, A_dist))}

    middle_pos = feature_pos_of_bond(od1, cg, N_dist)
    features |= {
        Feature('NEGC', above(middle_pos, cp, N_width)),
        Feature('NEGC', below(middle_pos, cp, N_width)),
        Feature('HDON', feature_pos_of_bond(od2, cg, A_dist)),
    }

    middle_pos = feature_pos_of_bond(od2, cg, N_dist)
    features |= {
        Feature('NEGC', above(middle_pos, cp, N_width)),
        Feature('NEGC', below(middle_pos, cp, N_width)),
        Feature('NEGC', feature_pos_of_bond(cg, cb, N_dist)),
    }

    return features


def features_from_sidechain_sulfur(residue):
    features = set()
    swapped_res = add_hydrogens2sulfur_as_carbon(residue)
    for c in sidechain_carbons(swapped_res):
        for hyd in bonded_hydrogens(c):
            feature_pos = feature_pos_of_bond(hyd, c, H_dist)
            feature = Feature('HDON', feature_pos)
            features.add(feature)
    return features


def features_from_cysteine_sidechain(residue):
    sg = residue.atom(name='SG')
    cb = residue.atom(name='CB')

    feature = Feature('HDON', feature_pos_of_bond(sg, cb, H_dist))
    features = {feature}

    features |= features_from_sidechain_sulfur(residue)

    return features


def features_from_histidines_sidechain(residue):
    ne2 = residue.atom(name='NE2')
    ce1 = residue.atom(name='CE1')
    cd2 = residue.atom(name='CD2')

    cp = center_of_triangle(ne2, cd2, ce1)

    features = set()
    ring_atom_names = {'NE2', 'CE1', 'CD2', 'ND1'}
    ring_atoms = [a for a in residue.atoms() if a.name() in ring_atom_names]
    ring_hydrogens = set()
    for ring_atom in ring_atoms:
        ring_hydrogens.update(bonded_hydrogens(ring_atom))
    if len(ring_hydrogens) == 4:
        for ring_atom in ring_atoms:
            hyd = bonded_hydrogens(ring_atom)[0]
            middle_pos = feature_pos_of_bond(hyd, ring_atom, P_dist)
            features |= {
                Feature('POSC', above(middle_pos, cp, P_width)),
                Feature('POSC', below(middle_pos, cp, P_width)),
            }

    ring_nitrogen_names = {'NE2', 'ND1'}
    ring_nitrogens = [a for a in residue.atoms() if a.name() in ring_nitrogen_names]
    for nitrogen in ring_nitrogens:
        hydrogens = bonded_hydrogens(nitrogen)
        for hydrogen in hydrogens:
            pos = feature_pos_of_bond(hydrogen, nitrogen, O_dist)
            feature = Feature('HACC', pos)
            features.add(feature)
        if not hydrogens:
            logging.warning("TODO: Not adding hydrogens to aromatic nitrogen, , less features will be generated until this is implemented")

    ring_center = center_of_ring(residue, {'CG', 'NE2', 'CE1', 'CD2', 'ND1'})
    features |= {
        Feature('AROM', above(ring_center, cp, Rp_width)),
        Feature('AROM', below(ring_center, cp, Rp_width)),
    }

    return features


def features_from_glutamicacids_sidechain(residue):
    cd = residue.atom(name='CD')
    oe1 = residue.atom(name='OE1')
    oe2 = residue.atom(name='OE2')

    features = set()
    cp = center_of_triangle(oe1, cd, oe2)

    feature = Feature('HDON', feature_pos_of_bond(oe1, cd, A_dist))
    features.add(feature)

    middle_pos = feature_pos_of_bond(oe1, cd, N_dist)
    features |= {
        Feature('NEGC', above(middle_pos, cp, N_width)),
        Feature('NEGC', below(middle_pos, cp, N_width)),
    }

    features.add(Feature('HDON', feature_pos_of_bond(oe2, cd, A_dist)))

    middle_pos = feature_pos_of_bond(oe2, cd, N_dist)
    features |= {
        Feature('NEGC', above(middle_pos, cp, N_width)),
        Feature('NEGC', below(middle_pos, cp, N_width)),
    }

    cg = residue.atom(name='CG')
    features.add(Feature('NEGC', feature_pos_of_bond(cd, cg, N_dist)))

    return features


def features_from_glutamines_sidechain(residue):
    cd = residue.atom(name='CD')
    oe1 = residue.atom(name='OE1')
    features = {Feature('HDON', feature_pos_of_bond(oe1, cd, A_dist))}

    ne2 = residue.atom(name='NE2')
    for hydrogen in bonded_hydrogens(ne2):
        feature = Feature('HACC', feature_pos_of_bond(hydrogen, ne2, O_dist))
        features.add(feature)

    return features


def features_from_isoleucine_sidechain(residue):
    features = set()

    cg1 = residue.atom(name='CG1')
    for hyd in bonded_hydrogens(cg1):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg1, H_dist)))

    cd1 = residue.atom(name='CD1')
    features.add(Feature('LIPO', feature_pos_of_bond(cd1, cg1, H_dist)))
    for hyd in bonded_hydrogens(cd1):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cd1, H_dist)))

    cg2 = residue.atom(name='CG2')
    for hyd in bonded_hydrogens(cg2):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg2, H_dist)))

    cb = residue.atom(name='CB')
    features.add(Feature('LIPO', feature_pos_of_bond(cg2, cb, H_dist)))
    for hyd in bonded_hydrogens(cb):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cb, H_dist)))

    return features


def features_from_leucine_sidechain(residue):
    features = set()

    cg = residue.atom(name='CG')
    for hyd in bonded_hydrogens(cg):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg, H_dist)))

    cd1 = residue.atom(name='CD1')
    features.add(Feature('LIPO', feature_pos_of_bond(cd1, cg, H_dist)))
    for hyd in bonded_hydrogens(cd1):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cd1, H_dist)))

    cb = residue.atom(name='CB')
    for hyd in bonded_hydrogens(cb):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cb, H_dist)))

    cd2 = residue.atom(name='CD2')
    features.add(Feature('LIPO', feature_pos_of_bond(cd2, cg, H_dist)))
    for hyd in bonded_hydrogens(cb):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cd2, H_dist)))

    return features


def features_from_lysine_donor_charged(residue):
    ce = residue.atom(name='CE')
    nz = residue.atom(name='NZ')

    features = {Feature('POSC', feature_pos_of_bond(nz, ce, P_dist))}

    for hyd in bonded_hydrogens(nz):
        features |= {
            Feature('HACC', feature_pos_of_bond(hyd, nz, O_dist)),
            Feature('POSC', feature_pos_of_bond(hyd, nz, P_dist)),
        }

    return features


def features_from_lysine_hydropobe(residue):
    features = set()

    cb = residue.atom(name='CB')
    for hyb in bonded_hydrogens(cb):
        features.add(Feature('LIPO', feature_pos_of_bond(hyb, cb, H_dist)))

    cg = residue.atom(name='CG')
    for hyb in bonded_hydrogens(cg):
        features.add(Feature('LIPO', feature_pos_of_bond(hyb, cg, H_dist)))

    cd = residue.atom(name='CD')
    for hyb in bonded_hydrogens(cd):
        features.add(Feature('LIPO', feature_pos_of_bond(hyb, cd, H_dist)))

    return features


def features_from_methionine_sidechain(residue):
    ce = residue.atom(name='CE')
    sd = residue.atom(name='SD')

    feature = Feature('HDON', feature_pos_of_bond(ce, sd, H_dist))
    features = {feature}
    features |= features_from_sidechain_sulfur(residue)
    return features


def features_from_phenylalanine_sidechain(residue):
    cg = residue.atom(name='CG')
    cd1 = residue.atom(name='CD1')
    cz = residue.atom(name='CZ')

    features = set()

    cp = center_of_triangle(cg, cd1, cz)
    ring_atoms = set()
    for atom in residue.atoms():
        if atom.name() in {'N', 'CA', 'C', 'O', 'CB'} or atom.element() == 'H':
            continue
        ring_atoms.add(atom)
    ring_center = (
        sum([a.location()[0] for a in ring_atoms]) / len(ring_atoms),
        sum([a.location()[1] for a in ring_atoms]) / len(ring_atoms),
        sum([a.location()[2] for a in ring_atoms]) / len(ring_atoms),
    )
    above_pos = (
        ring_center[0] + cp[0] * Rp_width,
        ring_center[1] + cp[1] * Rp_width,
        ring_center[2] + cp[2] * Rp_width,
    )
    features.add(Feature('AROM', above_pos))
    below_pos = (
        ring_center[0] - cp[0] * Rp_width,
        ring_center[1] - cp[1] * Rp_width,
        ring_center[2] - cp[2] * Rp_width,
    )
    features.add(Feature('AROM', below_pos))
    above_pos = (
        ring_center[0] + cp[0] * H_dist,
        ring_center[1] + cp[1] * H_dist,
        ring_center[2] + cp[2] * H_dist,
    )
    features.add(Feature('LIPO', above_pos))
    below_pos = (
        ring_center[0] - cp[0] * H_dist,
        ring_center[1] - cp[1] * H_dist,
        ring_center[2] - cp[2] * H_dist,
    )
    features.add(Feature('LIPO', below_pos))

    for ring_atom in ring_atoms:
        hyds = bonded_hydrogens(ring_atom)
        if hyds:
            hyd = hyds[0]
            features.add(Feature('AROM', feature_pos_of_bond(hyd, ring_atom, Rt_dist)))
            features.add(Feature('LIPO', feature_pos_of_bond(hyd, ring_atom, H_dist)))

    return features


def features_from_proline_sidechain(residue):
    features = set()
    carbons = residue.atoms(name='CD') | residue.atoms(name='CG') | residue.atoms(name='CB')
    for carbon in carbons:
        for hyd in bonded_hydrogens(carbon):
            features.add(Feature('LIPO', feature_pos_of_bond(hyd, carbon, H_dist)))
    return features


def features_from_hydroxyl_sidechain(residue, oxygen_name, hydrogen_name, carbon_name):
    deprotonated = len(residue.atoms(name=hydrogen_name)) == 0
    cb = residue.atom(name=carbon_name)
    og = residue.atom(name=oxygen_name)

    features = set()
    if not og:
        return features
    og_hyds = bonded_hydrogens(og)
    if len(og_hyds) == 0:
        return features
    h = og_hyds[0]

    if not deprotonated:
        features.add(Feature('HACC', feature_pos_of_bond(h, og, O_dist)))

    feature_pos = feature_pos_of_bond_rotated(h, og, A_dist, 120, cb)
    features.add(Feature('HDON', feature_pos))
    feature_pos = feature_pos_of_bond_rotated(h, og, A_dist, -120, cb)
    features.add(Feature('HDON', feature_pos))
    return features


def features_from_threonine_sidechain(residue: Residue):
    cg2 = residue.atom(name='CG2')
    cb = residue.atom(name='CB')

    feature_pos = feature_pos_of_bond(cg2, cb, H_dist)
    features = {Feature('LIPO', feature_pos)}

    for hyd in bonded_hydrogens(cg2):
        feature_pos = feature_pos_of_bond(hyd, cg2, H_dist)
        feature = Feature('LIPO', feature_pos)
        features.add(feature)

    return features


def features_from_sidechain_donor(residue: Residue):
    ne1 = residue.atom(name='NE1')
    he1 = [a for a in ne1.bonded_atoms() if a.element() == 'H'][0]

    feature_pos = feature_pos_of_bond(he1, ne1, O_dist)
    features = {Feature('HACC', feature_pos)}

    return features


def features_from_tryptophan_sidechain(residue: Residue):
    features = set()

    ne1 = residue.atom(name='NE1')
    cg = residue.atom(name='CG')
    ce2 = residue.atom(name='CE2')

    cp = center_of_triangle(ne1, cg, ce2)

    ring_center = center_of_ring(residue, {'CG', 'CD1', 'CD2', 'NE1', 'CE2'})
    above_pos = (
        ring_center[0] + cp[0] * Rp_width,
        ring_center[1] + cp[1] * Rp_width,
        ring_center[2] + cp[2] * Rp_width,
    )
    features.add(Feature('AROM', above_pos))
    below_pos = (
        ring_center[0] - cp[0] * Rp_width,
        ring_center[1] - cp[1] * Rp_width,
        ring_center[2] - cp[2] * Rp_width,
    )
    features.add(Feature('AROM', below_pos))

    ring_atom_names = {'CE3', 'CZ2', 'CZ3', 'CH2'}
    ring_center = center_of_ring(residue, ring_atom_names)
    above_pos = (
        ring_center[0] + cp[0] * Rp_width,
        ring_center[1] + cp[1] * Rp_width,
        ring_center[2] + cp[2] * Rp_width,
    )
    features.add(Feature('AROM', above_pos))
    below_pos = (
        ring_center[0] - cp[0] * Rp_width,
        ring_center[1] - cp[1] * Rp_width,
        ring_center[2] - cp[2] * Rp_width,
    )
    features.add(Feature('AROM', below_pos))

    above_pos = (
        ring_center[0] + cp[0] * H_dist,
        ring_center[1] + cp[1] * H_dist,
        ring_center[2] + cp[2] * H_dist,
    )
    features.add(Feature('LIPO', above_pos))
    below_pos = (
        ring_center[0] - cp[0] * H_dist,
        ring_center[1] - cp[1] * H_dist,
        ring_center[2] - cp[2] * H_dist,
    )
    features.add(Feature('LIPO', below_pos))

    for r in atoms_of_ring(residue, ring_atom_names):
        for h in bonded_hydrogens(r):
            features.add(Feature('AROM', feature_pos_of_bond(h, r, Rt_dist)))
            features.add(Feature('LIPO', feature_pos_of_bond(h, r, H_dist)))

    return features


def features_from_tyrosine_sidechain_ring(residue: Residue):
    features = set()

    cz = residue.atom(name='CZ')
    cg = residue.atom(name='CG')
    ce1 = residue.atom(name='CE1')

    if not (cz and cg and ce1):
        return features

    cp = center_of_triangle(cz, cg, ce1)

    ring_center = center_of_ring(residue, {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'})
    above_pos = (
        ring_center[0] + cp[0] * Rp_width,
        ring_center[1] + cp[1] * Rp_width,
        ring_center[2] + cp[2] * Rp_width,
    )
    features.add(Feature('AROM', above_pos))
    below_pos = (
        ring_center[0] - cp[0] * Rp_width,
        ring_center[1] - cp[1] * Rp_width,
        ring_center[2] - cp[2] * Rp_width,
    )
    features.add(Feature('AROM', below_pos))

    ring_center = center_of_ring(residue, {'CG', 'CD1', 'CD2'})
    above_pos = (
        ring_center[0] + cp[0] * H_dist,
        ring_center[1] + cp[1] * H_dist,
        ring_center[2] + cp[2] * H_dist,
    )
    features.add(Feature('LIPO', above_pos))
    below_pos = (
        ring_center[0] - cp[0] * H_dist,
        ring_center[1] - cp[1] * H_dist,
        ring_center[2] - cp[2] * H_dist,
    )
    features.add(Feature('LIPO', below_pos))

    ring_atom_names = {'CD1', 'CD2', 'CE1', 'CE2'}
    for r in atoms_of_ring(residue, ring_atom_names):
        for h in bonded_hydrogens(r):
            features |= {
                Feature('AROM', feature_pos_of_bond(h, r, Rt_dist)),
                Feature('LIPO', feature_pos_of_bond(h, r, H_dist)),
            }

    return features


def features_from_valine_sidechain(residue: Residue):
    features = set()

    cb = residue.atom(name='CB')
    cg1 = residue.atom(name='CG1')
    cg2 = residue.atom(name='CG2')

    for hyd in bonded_hydrogens(cb):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cb, H_dist)))

    features.add(Feature('LIPO', feature_pos_of_bond(cg1, cb, H_dist)))
    for hyd in bonded_hydrogens(cg1):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg1, H_dist)))

    features.add(Feature('LIPO', feature_pos_of_bond(cg2, cb, H_dist)))
    for hyd in bonded_hydrogens(cg1):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg2, H_dist)))

    return features


def features_from_alanine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_alanine_sidechain(residue),
    )


def features_from_arginine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_sidechain_hydrophobes(residue),
        features_from_sidechain_donors(residue),
        features_from_sidechain_positives(residue),
    )


def features_from_asparagine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_asparagine_sidechain(residue),
    )


def features_from_asparticacid(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_asparticacids_sidechain(residue),
    )


def features_from_cysteine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_cysteine_sidechain(residue),
    )


def features_from_histidine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_histidines_sidechain(residue),
    )


def features_from_glutamicacid(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_glutamicacids_sidechain(residue),
    )


def features_from_glutamine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_glutamines_sidechain(residue),
    )


def features_from_glycine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
    )


def features_from_isoleucine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_isoleucine_sidechain(residue),
    )


def features_from_leucine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_leucine_sidechain(residue),
    )


def features_from_lysine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_lysine_donor_charged(residue),
        features_from_lysine_hydropobe(residue),
    )


def features_from_methionine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_methionine_sidechain(residue),
    )


def features_from_phenylalanine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_phenylalanine_sidechain(residue),
    )


def features_from_proline(residue: Residue, _ligand):
    return set().union(
        features_from_backbone_carbonyl(residue),
        features_from_proline_sidechain(residue),
    )


def features_from_serine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_hydroxyl_sidechain(residue, 'OG', 'HG', 'CB'),
    )


def features_from_threonine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_hydroxyl_sidechain(residue, 'OG1', 'HG1', 'CB'),
        features_from_threonine_sidechain(residue),
    )


def features_from_tryptophan(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_sidechain_donor(residue),
        features_from_tryptophan_sidechain(residue),
    )


def features_from_tyrosine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_hydroxyl_sidechain(residue, 'OH', 'HH', 'CZ'),
        features_from_tyrosine_sidechain_ring(residue),
    )


def features_from_valine(residue: Residue, ligand):
    return set().union(
        features_from_backbone_amine(residue, ligand),
        features_from_backbone_carbonyl(residue),
        features_from_valine_sidechain(residue),
    )
