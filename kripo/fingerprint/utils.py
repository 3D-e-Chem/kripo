from typing import List, Generator, Tuple


def calc_bins(fp_width, fp_multiplier, max_bins, max_distance):
    """Calculates distance bins

    Args:
        fp_width (float):
        fp_multiplier (float):
        max_bins (int): Number of bins
        max_distance (float): Maximum distance

    Returns:
        List[float]: List of bin thresholds

    """
    bins = []
    value = float(0)
    previous_step = fp_width
    previous_value = fp_width
    bins.append(fp_width)
    while float(value) < max_distance and len(bins) < max_bins:
        step = previous_step * fp_multiplier
        value = previous_value + step
        bins.append(value)
        previous_value = value
        previous_step = step
    return bins


def bin_distance(distance, bins):
    """Bin the distance

    Args:
        distance (float): A distance
        bins (List[float]): The bins

    Returns:
        int: Index of bin in which distance is binned
    """
    nr_bins = len(bins)
    bin_id = nr_bins
    for i in range(nr_bins):
        if distance < bins[i]:
            bin_id = i
            break
    return bin_id


def calculate_distance_matrix(ordered_features) -> List[List[float]]:
    """Calculate distances between list of features

    Args:
        ordered_features (List[Feature]):

    Returns:
        Where row/column list indices are same a ordered_features indices and
            value is the distance between the row and column feature
    """
    n = len(ordered_features)
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            dist = ordered_features[i].distance(ordered_features[j])
            row.append(dist)
        matrix.append(row)

    return matrix


def fuzzy_offsets(factor: int=1, shape='all') -> Generator[Tuple[int, int, int], None, None]:
    """Generator for the fuzzy offsets

    Args:
        factor: The amount of fuzz to add.
        shape: Shape of offsets

            * all,
            * one,
            * v1,

    Yields:
        The offsets as x, y, z coordinates
    """
    if factor >= 0:
        for i in range((0 - factor), factor + 1):
            for j in range((0 - factor), factor + 1):
                for k in range((0 - factor), factor + 1):
                    yield i, j, k
    elif factor == -1:
        theset = [(1, 2, 1), (1, 1, -1), (1, -1, -1), (1, 0, -1), (0, -1, -1), (-1, -1, -1), (1, 2, 0)]
        for i, j, k in theset:
            yield i, j, k
    elif factor == -2:
        theset = [(0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
        for i, j, k in theset:
            yield i, j, k
    else:
        raise ValueError('Invalid fuzzy_factor')
