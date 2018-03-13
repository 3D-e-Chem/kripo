from atomium.files.pdb import Pdb

from kripo.pdb import ligands


def test_code(pdb_3heg: Pdb):
    assert pdb_3heg.code() == '3HEG'


def test_ligands(pdb_3heg: Pdb, ligand_expo_dict_fixture):
    ligand_names = [l.name() for l in ligands(pdb_3heg, ligand_expo_dict_fixture)]

    assert ligand_names == ['BAX']


def test_ligands__firstinstance(pdb_5is0: Pdb, ligand_expo_dict_fixture):
    ligand_names = [l.name() for l in ligands(pdb_5is0, ligand_expo_dict_fixture)]

    assert ligand_names == ['6ET']
