from atomium.structures.chains import Site

from kripo.site import chain_of_site


def test_chain_of_site(site_3heg_bax: Site):
    chain = chain_of_site(site_3heg_bax)

    assert chain == 'A'
