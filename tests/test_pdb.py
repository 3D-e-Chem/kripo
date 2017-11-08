from kripo.pdb import Pdb


def test_code(pdb_3heg: Pdb):
    assert pdb_3heg.code() == '3HEG'


def test_ligands(pdb_3heg: Pdb):
    ligands = pdb_3heg.ligands()
    ligand_names = [l.name() for l in ligands]

    assert ligand_names == ['BAX']


def test_ligands__firstinstance(pdb_5is0: Pdb):
    ligands = pdb_5is0.ligands()
    ligand_names = [l.name() for l in ligands]

    assert ligand_names == ['6ET']
