from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.structures import Model

from kripo.pdb import Pdb, pdb_from_file, ligands
from kripo.protonate import protonate_protein, protonate_ligand, protonate


def test_protonate_protein(orig_pdb_3heg: Pdb):
    block = orig_pdb_3heg.model().chain('A').to_file_string('pdb')
    # Verify input contains no hydrogens
    assert len(orig_pdb_3heg.model().chain('A').atoms(element='H')) == 0

    hblock, err = protonate_protein(block)

    hpdb = pdb_string2pdb(hblock)
    assert len(hpdb.model().chain('A').atoms(element='H')) == 2702


def test_protonate_ligand(orig_pdb_3heg: Pdb):
    # Create model which only contains ligand molecules
    ligands_model = Model()
    bax = orig_pdb_3heg.model().molecule(name='BAX')
    ligands_model.add_molecule(bax)
    block = ligands_model.to_file_string('pdb')
    # Verify input contains no hydrogens
    assert len(ligands_model.atoms(element='H')) == 0

    hblock = protonate_ligand(block)

    hpdb = pdb_string2pdb(hblock)
    hmodel = hpdb.model()
    assert len(hmodel.atoms(element='H')) == 16


def test_protonate_pdb_3heg(orig_pdb_3heg: Pdb):
    protonated = protonate(orig_pdb_3heg)
    assert len(protonated.model().atoms(element='H')) == 2936


def test_protonate_pdb_3rze():
    unprotonated = pdb_from_file('tests/fixtures/3RZE.pdb', hydrogenate=False, clean=False)
    protonated = protonate(unprotonated)
    assert len(protonated.model().atoms(element='H')) == 3655


def test_protonate_pdb_5is0():
    unprotonated = pdb_from_file('tests/fixtures/5IS0.pdb', hydrogenate=False, clean=False)
    protonated = protonate(unprotonated)
    assert len(protonated.model().atoms(element='H')) == 11072


def test_protonate_pdb_2muv():
    unprotonated = pdb_from_file('tests/fixtures/2MUV.pdb', hydrogenate=False, clean=False)
    protonated = protonate(unprotonated)
    assert len(protonated.model().atoms(element='H')) == 982


def pdb_string2pdb(block):
    return pdb_dict_to_pdb(pdb_string_to_pdb_dict(block))
