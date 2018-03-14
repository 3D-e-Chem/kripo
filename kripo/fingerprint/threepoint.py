import pkg_resources

from intbitset import intbitset

from .bitinfo import load_bitinfo
from .utils import calc_bins, bin_distance, calculate_distance_matrix, fuzzy_offsets

"""Constants and methods for 3 point pharmacophore fingerprints"""

"""List of distance bin thresholds"""
BINS = calc_bins(fp_width=0.8, fp_multiplier=1, max_bins=20, max_distance=float(23))

"""Dictionary of that maps a 3 point pharmacophore to a fingerprint bit position

The key is a string formatted as `XYZabc`. Where `XYZ` are the pharmacophore type short names in alphabetical order
and where `abc` are the binned charified distances between the features in shortest distance first order.
The shortest bin distance gets char `a` and each next bin is the next char in the alphabet.
"""
BIT_INFO = load_bitinfo(pkg_resources.resource_filename('kripo', 'data/PHARMACKEY_3PFP_25BINS_6FEATURES_new.txt.bz2'))

"""Pharmacophore feature type name to fingerprint bit short name"""
FEATURE2BIT = {
    'LIPO': 'H',
    'POSC': 'P',
    'NEGC': 'N',
    'HDON': 'O',
    'HACC': 'A',
    'AROM': 'R'
}


def from_pharmacophore(pharmacophore, subs=True, fuzzy_factor=1, fuzzy_shape='all') -> intbitset:
    """Build a fingerprint from a pharmacophore

    Args:
        pharmacophore (Pharmacophore): The pharmacophore
        subs (bool): Include bits for 1 point and 2 points
        fuzzy_factor (int): Number of bins below/above actual bin to include in fingerprint
        fuzzy_shape (str): Shape of fuzzying

    Returns:
        Fingerprint
    """
    ordered_features = sorted(list(pharmacophore.features), key=lambda f: FEATURE2BIT[f.kind])
    nr_features = len(ordered_features)
    nr_bins = len(BINS)
    nr_bins_three_point = nr_bins
    if fuzzy_shape == 'v1':
        # v1 offsets fail bounds check, so bin indices outside nr_bins are allowed
        nr_bins_three_point += fuzzy_factor * 2

    dist_matrix = calculate_distance_matrix(ordered_features)

    bits = set()

    if subs:
        for a in ordered_features:
            bit_info = FEATURE2BIT[a.kind]
            bit_index = BIT_INFO[bit_info]
            bits.add(bit_index)

    if subs:
        pair_fuzzy_factor = fuzzy_factor
        if pair_fuzzy_factor < 0:
            pair_fuzzy_factor = 1
        for a in range(nr_features - 1):
            for b in range(a + 1, nr_features):
                distance = dist_matrix[a][b]
                bin_id = bin_distance(distance, BINS)
                ids = [
                    FEATURE2BIT[ordered_features[a].kind],
                    FEATURE2BIT[ordered_features[b].kind]
                ]
                bit_info = ids[0] + chr(bin_id + 97) + ids[1]
                bit_index = BIT_INFO[bit_info]
                bits.add(bit_index)

                for i in range((0 - pair_fuzzy_factor), pair_fuzzy_factor + 1):
                    if nr_bins < bin_id + i or bin_id + i < 0:
                        continue

                    bit_info = ids[0] + chr(bin_id + 97 + i) + ids[1]
                    bit_index = BIT_INFO[bit_info]
                    bits.add(bit_index)

    offsets = tuple(set(fuzzy_offsets(fuzzy_factor, fuzzy_shape)))
    for a in range(nr_features - 2):
        for b in range(a + 1, nr_features - 1):
            for c in range(b + 1, nr_features):
                distances = [(
                    bin_distance(dist_matrix[a][b], BINS),
                    FEATURE2BIT[ordered_features[c].kind],
                ), (
                    bin_distance(dist_matrix[a][c], BINS),
                    FEATURE2BIT[ordered_features[b].kind],
                ), (
                    bin_distance(dist_matrix[b][c], BINS),
                    FEATURE2BIT[ordered_features[a].kind],
                )]
                distances.sort()

                bit_info = distances[0][1] + distances[1][1] + distances[2][1] \
                    + chr(distances[0][0] + 97) \
                    + chr(distances[1][0] + 97) \
                    + chr(distances[2][0] + 97)
                bit_index = BIT_INFO[bit_info]
                bits.add(bit_index)

                for i, j, k in offsets:
                    # test if is bit outside bins
                    bin_i = distances[0][0] + i
                    bin_j = distances[1][0] + j
                    bin_k = distances[2][0] + k

                    if bin_i < 0 or bin_i >= nr_bins_three_point or \
                            bin_j < 0 or bin_j >= nr_bins_three_point or \
                            bin_k < 0 or bin_k >= nr_bins_three_point:
                        continue

                    fdistances = [(
                        bin_i,
                        distances[0][1],
                    ), (
                        bin_j,
                        distances[1][1],
                    ), (
                        bin_k,
                        distances[2][1],
                    )]

                    fdistances.sort()

                    bit_info = fdistances[0][1] + fdistances[1][1] + fdistances[2][1] \
                        + chr(fdistances[0][0] + 97) \
                        + chr(fdistances[1][0] + 97) \
                        + chr(fdistances[2][0] + 97)
                    bit_index = BIT_INFO[bit_info]
                    bits.add(bit_index)

    return intbitset(bits)
