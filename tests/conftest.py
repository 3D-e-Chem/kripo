import pytest

from kripo.ligand import Ligand
from kripo.pdb import pdb_from_file, ligands


def prep_pdb(opath, npath):
    """Speed up tests by using already hydrogenated and cleaned pdbs"""
    pdb = pdb_from_file(opath)
    pdb.save(npath)


@pytest.fixture(scope="module")
def pdb_3heg():
    filename = 'tests/fixtures/3HEG.prepped.pdb'
    return pdb_from_file(filename, hydrogenate=False, clean=False)


@pytest.fixture(scope="module")
def pdb_5is0():
    # electron microscopy
    filename = 'tests/fixtures/5IS0.prepped.pdb'
    return pdb_from_file(filename, hydrogenate=False, clean=False)


@pytest.fixture(scope="module")
def pdb_3rze():
    # multiple ligands
    filename = 'tests/fixtures/3RZE.prepped.pdb'
    return pdb_from_file(filename, hydrogenate=False, clean=False)


@pytest.fixture(scope="module")
def pdb_2muv():
    # NMR and multi model
    filename = 'tests/fixtures/2MUV.prepped.pdb'
    return pdb_from_file(filename, hydrogenate=False, clean=False)


@pytest.fixture(scope="module")
def ligand_3heg_bax(pdb_3heg) -> Ligand:
    return ligands(pdb_3heg)[0]


@pytest.fixture(scope="module")
def site_3heg_bax(ligand_3heg_bax: Ligand):
    return ligand_3heg_bax.site()


@pytest.fixture(scope="module")
def fragments_3heg_bax(ligand_3heg_bax: Ligand):
    return ligand_3heg_bax.fragments()


@pytest.fixture(scope="module")
def fragment1_3heg_bax(fragments_3heg_bax):
    return fragments_3heg_bax[0]


@pytest.fixture(scope="module")
def fragment2_3heg_bax(fragments_3heg_bax):
    return fragments_3heg_bax[1]


@pytest.fixture(scope="module")
def yasara_fragment2_3heg_bax():
    filename = 'tests/fixtures/3HEG.frag2.pdb'
    pdb = pdb_from_file(filename)
    ligand = ligands(pdb)[0]
    frags = ligand.fragments()
    return frags[0]
