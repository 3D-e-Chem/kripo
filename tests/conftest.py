import pytest

from kripo.ligand import Ligand
from kripo.pdb import Pdb


@pytest.fixture(scope="module")
def pdb_3heg():
    filename = 'tests/fixtures/3HEG.pdb'
    return Pdb(filename)


@pytest.fixture(scope="module")
def pdb_5is0():
    # electron microscopy
    filename = 'tests/fixtures/5IS0.pdb'
    return Pdb(filename)


@pytest.fixture(scope="module")
def pdb_3rze():
    # multiple ligands
    filename = 'tests/fixtures/3RZE.pdb'
    return Pdb(filename)


@pytest.fixture(scope="module")
def pdb_2muv():
    # NMR and multi model
    filename = 'tests/fixtures/2MUV.pdb'
    return Pdb(filename)


@pytest.fixture(scope="module")
def ligand_3heg_bax(pdb_3heg: Pdb):
    return list(pdb_3heg.ligands())[0]


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
