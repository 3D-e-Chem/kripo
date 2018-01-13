import logging

from atomium.structures import Residue

from .feature import Feature
from .utils import feature_pos_of_bond, feature_pos_of_bond_rotated, bonded_hydrogens, atoms_by_name, \
    center_of_atoms_by_name, sidechain_nitrogens, sidechain_carbons, add_hydrogens2sulfur_as_carbon
from .vector import center_of_triangle, above, below

H_dist = 0.8  # HYDROPHOBE
Rp_width = 4.0  # PLACEMENT OF PI AROMATIC FEATURES
Rt_dist = -1  # ATOMATIC T STACKING
A_dist = 0.8  # ACCEPTORS
O_dist = 0.2  # DONORS
P_dist = 0.3  # POSITIVE CHARGES
P_width = 0.0  # POSITIVE CHARGES
N_dist = 1.4  # NEGATIVE CHARGES
N_width = 1.0  # NEGATIVE CHARGES


def features_from_backbone_amine(residue: Residue):
    if O_dist < 0:
        return set()
    nitrogen = residue.atom(name='N')
    try:
        hydrogen = bonded_hydrogens(nitrogen)[0]
        feature_pos = feature_pos_of_bond(hydrogen, nitrogen, O_dist)
        return {Feature('HDON', feature_pos)}
    except IndexError:
        logging.warning('Skipping amine backbone feature as residue {0} has no N-H bond'.format(residue))
        return set()


def features_from_backbone_carbonyl(residue):
    if A_dist < 0:
        return set()
    oxygen = residue.atom(name='O')
    carbon = residue.atom(name='C')

    feature_pos = feature_pos_of_bond(oxygen, carbon, A_dist)
    feature = Feature('HACC', feature_pos)

    return {feature}


def features_from_alanine_sidechain(residue):
    if H_dist < 0:
        return set()
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
    features = set()
    if H_dist < 0:
        return features

    cb = residue.atom(name='CB')
    if cb:
        for hyd in bonded_hydrogens(cb):
            feature_pos = feature_pos_of_bond(hyd, cb, H_dist)
            feature = Feature('LIPO', feature_pos)
            features.add(feature)
    cg = residue.atom(name='CG')
    if cg:
        for hyd in bonded_hydrogens(cg):
            feature_pos = feature_pos_of_bond(hyd, cg, H_dist)
            feature = Feature('LIPO', feature_pos)
            features.add(feature)
    return features


def features_from_sidechain_donors(residue):
    features = set()
    if O_dist < 0:
        return features
    for nitrogen in sidechain_nitrogens(residue):
        for hydrogen in bonded_hydrogens(nitrogen):
            feature_pos = feature_pos_of_bond(hydrogen, nitrogen, O_dist)
            feature = Feature('HDON', feature_pos)
            features.add(feature)
    return features


def features_from_sidechain_positives(residue):
    features = set()
    if P_dist < 0:
        return features

    cz = residue.atom(name='CZ')
    nh1 = residue.atom(name='NH1')
    nh2 = residue.atom(name='NH2')

    if cz and nh1 and nh2:
        cp = center_of_triangle(nh1, cz, nh2)

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

    features = set()
    if A_dist >= 0 and od1 and cg:
        feature_pos = feature_pos_of_bond(od1, cg, A_dist)
        feature = Feature('HACC', feature_pos)
        features.add(feature)

    if O_dist >= 0 and nd2:
        for hyd in bonded_hydrogens(nd2):
            feature_pos = feature_pos_of_bond(hyd, nd2, O_dist)
            feature = Feature('HDON', feature_pos)
            features.add(feature)

    return features


def features_from_asparticacids_sidechain(residue):
    cb = residue.atom(name='CB')
    cg = residue.atom(name='CG')
    od1 = residue.atom(name='OD1')
    od2 = residue.atom(name='OD2')

    cp = center_of_triangle(od1, cg, od2)

    features = set()
    if A_dist >= 0 and cg:
        if od1:
            features.add(Feature('HACC', feature_pos_of_bond(od1, cg, A_dist)))
        if od2:
            features.add(Feature('HACC', feature_pos_of_bond(od2, cg, A_dist)))

    if N_dist >= 0 and cg:
        if od1:
            middle_pos = feature_pos_of_bond(od1, cg, N_dist)
            features |= {
                Feature('NEGC', above(middle_pos, cp, N_width)),
                Feature('NEGC', below(middle_pos, cp, N_width)),
            }

        if od2:
            middle_pos = feature_pos_of_bond(od2, cg, N_dist)
            features |= {
                Feature('NEGC', above(middle_pos, cp, N_width)),
                Feature('NEGC', below(middle_pos, cp, N_width)),
                Feature('NEGC', feature_pos_of_bond(cg, cb, N_dist)),
            }

    return features


def features_from_sidechain_sulfur(residue):
    features = set()
    if H_dist < 0:
        return features
    for c in sidechain_carbons(residue):
        for hyd in bonded_hydrogens(c):
            feature_pos = feature_pos_of_bond(hyd, c, H_dist)
            feature = Feature('LIPO', feature_pos)
            features.add(feature)
    try:
        swapped_res = add_hydrogens2sulfur_as_carbon(residue)
        for s in swapped_res.atoms(element='S'):
            for hyd in bonded_hydrogens(s):
                feature_pos = feature_pos_of_bond(hyd, s, H_dist)
                feature = Feature('LIPO', feature_pos)
                features.add(feature)
    except ValueError:
        msg = 'Unable to determine plane for sulfur LIPO features of residue {0}'.format(residue)
        logging.warning(msg)
    return features


def features_from_cysteine_sidechain(residue):
    return features_from_sidechain_sulfur(residue)


def features_from_histidines_sidechain(residue):
    features = set()

    ring_nitrogen_names = {'NE2', 'ND1'}
    ring_nitrogens = [a for a in residue.atoms() if a.name() in ring_nitrogen_names]
    if O_dist >= 0:
        for nitrogen in ring_nitrogens:
            hydrogens = bonded_hydrogens(nitrogen)
            for hydrogen in hydrogens:
                pos = feature_pos_of_bond(hydrogen, nitrogen, O_dist)
                feature = Feature('HDON', pos)
                features.add(feature)
            if not hydrogens:
                logging.warning(
                    'TODO: Not adding hydrogens to aromatic nitrogen of HIS, ' +
                    'less features will be generated until this is implemented, ' +
                    'of residue {0}'.format(residue))

    ne2 = residue.atom(name='NE2')
    ce1 = residue.atom(name='CE1')
    cd2 = residue.atom(name='CD2')

    if not(ne2 and ce1 and cd2):
        return features

    cp = center_of_triangle(ne2, cd2, ce1)

    ring_atom_names = {'NE2', 'CE1', 'CD2', 'ND1'}
    ring_atoms = [a for a in residue.atoms() if a.name() in ring_atom_names]
    ring_hydrogens = set()
    for ring_atom in ring_atoms:
        ring_hydrogens.update(bonded_hydrogens(ring_atom))
    if len(ring_hydrogens) == 4 and P_width >= 0:
        for ring_atom in ring_atoms:
            hyd = bonded_hydrogens(ring_atom)[0]
            middle_pos = feature_pos_of_bond(hyd, ring_atom, P_dist)
            features |= {
                Feature('POSC', above(middle_pos, cp, P_width)),
                Feature('POSC', below(middle_pos, cp, P_width)),
            }

    if Rp_width >= 0:
        ring_center = center_of_atoms_by_name(residue, {'CG', 'NE2', 'CE1', 'CD2', 'ND1'})
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
    if A_dist >= 0:
        if oe1 and cd:
            features.add(Feature('HACC', feature_pos_of_bond(oe1, cd, A_dist)))
        if oe2 and cd:
            features.add(Feature('HACC', feature_pos_of_bond(oe2, cd, A_dist)))

    if N_dist >= 0 and oe1 and cd and oe2:
        cp = center_of_triangle(oe1, cd, oe2)
        middle_pos = feature_pos_of_bond(oe1, cd, N_dist)
        features |= {
            Feature('NEGC', above(middle_pos, cp, N_width)),
            Feature('NEGC', below(middle_pos, cp, N_width)),
        }

        middle_pos = feature_pos_of_bond(oe2, cd, N_dist)
        features |= {
            Feature('NEGC', above(middle_pos, cp, N_width)),
            Feature('NEGC', below(middle_pos, cp, N_width)),
        }

    if N_dist >= 0:
        cg = residue.atom(name='CG')
        if cd and cg:
            features.add(Feature('NEGC', feature_pos_of_bond(cd, cg, N_dist)))

    return features


def features_from_glutamines_sidechain(residue):
    cd = residue.atom(name='CD')
    oe1 = residue.atom(name='OE1')
    features = set()
    if A_dist >= 0 and cd and oe1:
        features.add(Feature('HACC', feature_pos_of_bond(oe1, cd, A_dist)))

    ne2 = residue.atom(name='NE2')
    if O_dist >= 0 and ne2:
        for hydrogen in bonded_hydrogens(ne2):
            feature = Feature('HDON', feature_pos_of_bond(hydrogen, ne2, O_dist))
            features.add(feature)

    return features


def features_from_isoleucine_sidechain(residue):
    features = set()
    if H_dist < 0:
        return features

    cg1 = residue.atom(name='CG1')
    if cg1:
        for hyd in bonded_hydrogens(cg1):
            features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg1, H_dist)))

    cd1 = residue.atom(name='CD1')
    if cg1 and cd1:
        features.add(Feature('LIPO', feature_pos_of_bond(cd1, cg1, H_dist)))
    if cd1:
        for hyd in bonded_hydrogens(cd1):
            # TODO check if LIPO on CD1 hydrogens should be added always, even if H_dist < 0
            features.add(Feature('LIPO', feature_pos_of_bond(hyd, cd1, H_dist)))

    cg2 = residue.atom(name='CG2')
    if cg2:
        for hyd in bonded_hydrogens(cg2):
            features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg2, H_dist)))

    cb = residue.atom(name='CB')
    if cb and cg2:
        features.add(Feature('LIPO', feature_pos_of_bond(cg2, cb, H_dist)))
    if cb:
        for hyd in bonded_hydrogens(cb):
            features.add(Feature('LIPO', feature_pos_of_bond(hyd, cb, H_dist)))

    return features


def features_from_leucine_sidechain(residue):
    features = set()
    if H_dist < 0:
        return features

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
    for hyd in bonded_hydrogens(cd2):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cd2, H_dist)))

    return features


def features_from_lysine_donor_charged(residue):
    ce = residue.atom(name='CE')
    nz = residue.atom(name='NZ')

    features = set()
    if P_dist >= 0 and ce and nz:
        features.add(Feature('POSC', feature_pos_of_bond(nz, ce, P_dist)))

    if nz:
        for hyd in bonded_hydrogens(nz):
            if O_dist >= 0:
                features.add(Feature('HDON', feature_pos_of_bond(hyd, nz, O_dist)))
            if P_dist >= 0:
                features.add(Feature('POSC', feature_pos_of_bond(hyd, nz, P_dist)))

    return features


def features_from_lysine_hydropobe(residue):
    features = set()
    if H_dist < 0:
        return features

    cb = residue.atom(name='CB')
    for hyb in bonded_hydrogens(cb):
        features.add(Feature('LIPO', feature_pos_of_bond(hyb, cb, H_dist)))

    cg = residue.atom(name='CG')
    if cg:
        for hyb in bonded_hydrogens(cg):
            features.add(Feature('LIPO', feature_pos_of_bond(hyb, cg, H_dist)))

    cd = residue.atom(name='CD')
    if cd:
        for hyb in bonded_hydrogens(cd):
            features.add(Feature('LIPO', feature_pos_of_bond(hyb, cd, H_dist)))

    return features


def features_from_methionine_sidechain(residue):
    features = features_from_sidechain_sulfur(residue)
    if H_dist >= 0:
        sg = residue.atom(name='SD')
        ce = residue.atom(name='CE')
        features.add(Feature('LIPO', feature_pos_of_bond(ce, sg, H_dist)))
    return features


def features_from_phenylalanine_sidechain(residue):
    cg = residue.atom(name='CG')
    cd1 = residue.atom(name='CD1')
    cz = residue.atom(name='CZ')

    features = set()

    cp = center_of_triangle(cg, cd1, cz)
    ring_atom_names = {'CG', 'CE2', 'CE1', 'CD2', 'CD1', 'CZ'}
    ring_center = center_of_atoms_by_name(residue, ring_atom_names)

    if Rp_width >= 0:
        features |= {
            Feature('AROM', above(ring_center, cp, Rp_width)),
            Feature('AROM', below(ring_center, cp, Rp_width)),
        }

    if H_dist >= 0:
        features |= {
            Feature('LIPO', above(ring_center, cp, H_dist)),
            Feature('LIPO', below(ring_center, cp, H_dist)),
        }

    ring_atoms = [a for a in residue.atoms() if a.name() in ring_atom_names]
    for ring_atom in ring_atoms:
        hyds = bonded_hydrogens(ring_atom)
        if hyds:
            hyd = hyds[0]
            if Rt_dist >= 0:
                features.add(Feature('AROM', feature_pos_of_bond(hyd, ring_atom, Rt_dist)))
            if H_dist >= 0:
                features.add(Feature('LIPO', feature_pos_of_bond(hyd, ring_atom, H_dist)))

    return features


def features_from_proline_sidechain(residue):
    features = set()
    if H_dist < 0:
        return features
    carbons = residue.atoms(name='CD') | residue.atoms(name='CG') | residue.atoms(name='CB')
    for carbon in carbons:
        for hyd in bonded_hydrogens(carbon):
            features.add(Feature('LIPO', feature_pos_of_bond(hyd, carbon, H_dist)))
    return features


def features_from_hydroxyl_sidechain(residue, oxygen_name, hydrogen_name, carbon_name):
    protonated = len(residue.atoms(name=hydrogen_name)) > 0
    cb = residue.atom(name=carbon_name)
    og = residue.atom(name=oxygen_name)

    features = set()
    if not og:
        logging.warning('Skipping hydroxyl sidechain features as residue {0} has no oxygen'.format(residue))
        return features
    og_hyds = bonded_hydrogens(og)
    if len(og_hyds) == 0:
        logging.warning('Skipping hydroxyl sidechain features as residue {0} has no O-H bond'.format(residue))
        return features
    h = og_hyds[0]

    if protonated and O_dist >= 0:
        features.add(Feature('HDON', feature_pos_of_bond(h, og, O_dist)))

    if A_dist >= 0:
        feature_pos = feature_pos_of_bond_rotated(h, og, A_dist, 120, cb)
        features.add(Feature('HACC', feature_pos))
        feature_pos = feature_pos_of_bond_rotated(h, og, A_dist, -120, cb)
        features.add(Feature('HACC', feature_pos))
    return features


def features_from_threonine_sidechain(residue: Residue):
    features = set()
    if H_dist < 0:
        return features

    cg2 = residue.atom(name='CG2')
    cb = residue.atom(name='CB')

    if cg2 and cb:
        feature_pos = feature_pos_of_bond(cg2, cb, H_dist)
        features = {Feature('LIPO', feature_pos)}

    if cg2:
        for hyd in bonded_hydrogens(cg2):
            feature_pos = feature_pos_of_bond(hyd, cg2, H_dist)
            feature = Feature('LIPO', feature_pos)
            features.add(feature)

    return features


def features_from_sidechain_donor(residue: Residue):
    if O_dist < 0:
        return set()
    ne1 = residue.atom(name='NE1')
    he1 = bonded_hydrogens(ne1)[0]

    feature_pos = feature_pos_of_bond(he1, ne1, O_dist)
    features = {Feature('HDON', feature_pos)}

    return features


def features_from_tryptophan_sidechain(residue: Residue):
    features = set()

    ne1 = residue.atom(name='NE1')
    cg = residue.atom(name='CG')
    ce2 = residue.atom(name='CE2')

    cp = center_of_triangle(ne1, cg, ce2)

    ring_center = center_of_atoms_by_name(residue, {'CG', 'CD1', 'CD2', 'NE1', 'CE2'})
    if Rp_width >= 0:
        features |= {
            Feature('AROM', above(ring_center, cp, Rp_width)),
            Feature('AROM', below(ring_center, cp, Rp_width)),
        }

    ring_atom_names = {'CE3', 'CZ2', 'CZ3', 'CH2'}
    ring_center = center_of_atoms_by_name(residue, ring_atom_names)
    if Rp_width >= 0:
        features |= {
            Feature('AROM', above(ring_center, cp, Rp_width)),
            Feature('AROM', below(ring_center, cp, Rp_width)),
        }
    if H_dist >= 0:
        features |= {
            Feature('LIPO', above(ring_center, cp, H_dist)),
            Feature('LIPO', below(ring_center, cp, H_dist)),
        }

    for r in atoms_by_name(residue, ring_atom_names):
        for h in bonded_hydrogens(r):
            if Rt_dist >= 0:
                features.add(Feature('AROM', feature_pos_of_bond(h, r, Rt_dist)))
            if H_dist >= 0:
                features.add(Feature('LIPO', feature_pos_of_bond(h, r, H_dist)))

    return features


def features_from_tyrosine_sidechain_ring(residue: Residue):
    features = set()

    cz = residue.atom(name='CZ')
    cg = residue.atom(name='CG')
    ce1 = residue.atom(name='CE1')

    if not (cz and cg and ce1):
        logging.warning('Skipping ring side chain features as residue {0} has no ring'.format(residue))
        return features

    cp = center_of_triangle(cz, cg, ce1)

    if Rp_width >= 0:
        ring_center = center_of_atoms_by_name(residue, {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'})
        features |= {
            Feature('AROM', above(ring_center, cp, Rp_width)),
            Feature('AROM', below(ring_center, cp, Rp_width)),
        }

    # TODO check if this is correct, it is weird to use negative number for position calc
    if H_dist < 0:
        ring_center = center_of_atoms_by_name(residue, {'CG', 'CD1', 'CD2'})
        features |= {
            Feature('LIPO', above(ring_center, cp, H_dist)),
            Feature('LIPO', below(ring_center, cp, H_dist)),
        }

    ring_atom_names = {'CD1', 'CD2', 'CE1', 'CE2'}
    for r in atoms_by_name(residue, ring_atom_names):
        for h in bonded_hydrogens(r):
            if Rt_dist >= 0:
                features.add(Feature('AROM', feature_pos_of_bond(h, r, Rt_dist)))
                if H_dist >= 0:
                    # TODO check if Rt_dist and H_dist should be >=0 or if H_dist >= 0 is enough to add feature
                    features.add(Feature('LIPO', feature_pos_of_bond(h, r, H_dist)))

    return features


def features_from_valine_sidechain(residue: Residue):
    features = set()
    if H_dist < 0:
        return features

    cb = residue.atom(name='CB')
    cg1 = residue.atom(name='CG1')
    cg2 = residue.atom(name='CG2')

    for hyd in bonded_hydrogens(cb):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cb, H_dist)))

    features.add(Feature('LIPO', feature_pos_of_bond(cg1, cb, H_dist)))
    for hyd in bonded_hydrogens(cg1):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg1, H_dist)))

    features.add(Feature('LIPO', feature_pos_of_bond(cg2, cb, H_dist)))
    for hyd in bonded_hydrogens(cg2):
        features.add(Feature('LIPO', feature_pos_of_bond(hyd, cg2, H_dist)))

    return features


def features_from_alanine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_alanine_sidechain(residue),
    )


def features_from_arginine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_sidechain_hydrophobes(residue),
        features_from_sidechain_donors(residue),
        features_from_sidechain_positives(residue),
    )


def features_from_asparagine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_asparagine_sidechain(residue),
    )


def features_from_asparticacid(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_asparticacids_sidechain(residue),
    )


def features_from_cysteine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_cysteine_sidechain(residue),
    )


def features_from_histidine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_histidines_sidechain(residue),
    )


def features_from_glutamicacid(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_glutamicacids_sidechain(residue),
    )


def features_from_glutamine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_glutamines_sidechain(residue),
    )


def features_from_glycine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
    )


def features_from_isoleucine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_isoleucine_sidechain(residue),
    )


def features_from_leucine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_leucine_sidechain(residue),
    )


def features_from_lysine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_lysine_donor_charged(residue),
        features_from_lysine_hydropobe(residue),
    )


def features_from_methionine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_methionine_sidechain(residue),
    )


def features_from_phenylalanine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_phenylalanine_sidechain(residue),
    )


def features_from_proline(residue: Residue):
    return set().union(
        features_from_backbone_carbonyl(residue),
        features_from_proline_sidechain(residue),
    )


def features_from_serine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_hydroxyl_sidechain(residue, 'OG', 'HG', 'CB'),
    )


def features_from_threonine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_hydroxyl_sidechain(residue, 'OG1', 'HG1', 'CB'),
        features_from_threonine_sidechain(residue),
    )


def features_from_tryptophan(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_sidechain_donor(residue),
        features_from_tryptophan_sidechain(residue),
    )


def features_from_tyrosine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_hydroxyl_sidechain(residue, 'OH', 'HH', 'CZ'),
        features_from_tyrosine_sidechain_ring(residue),
    )


def features_from_valine(residue: Residue):
    return set().union(
        features_from_backbone_amine(residue),
        features_from_backbone_carbonyl(residue),
        features_from_valine_sidechain(residue),
    )
