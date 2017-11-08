import bz2


def load_bitinfo(filename):
    """Loads a bit info file into a dictionary.

    Args:
        filename (str): Filename of bit info file

    Returns:
        Dict: Where key is info string and value is bit index.

    """
    data = {}
    with bz2.open(filename, 'rt') as f:
        for index, line in enumerate(f, 1):
            info = line.strip()
            data[info] = index
    return data
