from collections import Counter

from atomium.structures.chains import Site


def chain_of_site(site: Site) -> str:
    """Chain identifier of chain most site residues belong to

    A site can be in different chains.

    Args:
        site:

    Returns:

    """
    most_common = Counter([r.chain().chain_id() for r in site.residues()]).most_common(1)
    return most_common[0][0]
