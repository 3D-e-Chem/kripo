
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
