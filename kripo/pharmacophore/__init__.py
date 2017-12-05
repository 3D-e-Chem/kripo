from typing import Set

from atomium.structures.chains import Site

from .from_residues import features_from_alanine, features_from_arginine, \
    features_from_asparagine, features_from_asparticacid, features_from_cysteine, features_from_histidine, \
    features_from_glutamicacid, features_from_glutamine, features_from_glycine, features_from_isoleucine, \
    features_from_leucine, features_from_lysine, features_from_methionine, features_from_phenylalanine, \
    features_from_proline, features_from_serine, features_from_threonine, features_from_tryptophan, \
    features_from_tyrosine, features_from_valine
from .feature import Feature
from ..fingerprint.threepoint import from_pharmacophore
from ..fragment import Fragment


class NoFeatures(ValueError):
    pass


class Pharmacophore:
    """Pharmacophore

    Attributes:
        List[Feature]: List of features
    """
    def __init__(self, features):
        self.features = features

    def fingerprint(self):
        """Fingerprint of pharmacophore

        Returns:
            intbitset.intbitset: Fingerprint
        """
        return from_pharmacophore(self)

    def __repr__(self):
        return 'Pharmacophore({0})'.format(repr(self.features))


def from_fragment(fragment: Fragment):
    """Build a pharmacophore from a fragment

    Args:
        fragment (Fragment): the fragment

    Returns:
        Pharmacophore: the pharmacophore

    Raises:
        NoFeatures: When constructed pharmacophore has no features
    """
    features = from_site(fragment.site())

    features = filter_contact_features(fragment, features)

    if features:
        return Pharmacophore(features)
    else:
        raise NoFeatures


def from_site(site: Site):
    """Build pharmacophore features from a site

    Args:
        site (Site): The site

    Returns:
        List[Feature]: List of features

    """
    ligand = site.ligand()
    features = set()
    mappers = {
        'ALA': features_from_alanine,
        'ARG': features_from_arginine,
        'ASN': features_from_asparagine,
        'ASP': features_from_asparticacid,
        'CYS': features_from_cysteine,
        'HIS': features_from_histidine,
        'GLU': features_from_glutamicacid,
        'GLN': features_from_glutamine,
        'GLY': features_from_glycine,
        'ILE': features_from_isoleucine,
        'LEU': features_from_leucine,
        'LYS': features_from_lysine,
        'MET': features_from_methionine,
        'PHE': features_from_phenylalanine,
        'PRO': features_from_proline,
        'SER': features_from_serine,
        'THR': features_from_threonine,
        'TRP': features_from_tryptophan,
        'TYR': features_from_tyrosine,
        'VAL': features_from_valine,
    }
    for residue in site.residues():
        features |= mappers[residue.name()](residue)

    features = annihilate_neighbouring_donors_and_acceptors(features)

    return features


def annihilate_neighbouring_donors_and_acceptors(features: Set[Feature]):
    """Removes hydrogen donor/acceptor pairs which are 1.1 angstroms apart.

    Args:
        features (Set[Feature]): The features

    Returns:
        Set[Feature]: The features with possibly removed hydrogen donors and acceptors.

    """
    donors = [f for f in features if f.kind == 'HDON']
    acceptors = [f for f in features if f.kind == 'HACC']
    neighbours = set()
    for donor in donors:
        for acceptor in acceptors:
            if donor.distance(acceptor) < 1.1:
                neighbours.add(donor)
                neighbours.add(acceptor)

    return {f for f in features if f not in neighbours}


def filter_contact_features(fragment, features):
    """Filter on features that are in contact with fragment

    Args:
        fragment (Fragment): Fragment of ligand
        features (List[Feature):  List of features

    Returns:
        Set[Feature]: Features which are in contact with fragment
    """
    return {f for f in features if f.in_contact_with(fragment)}
