from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from rdkit.Chem import Mol

from kripo.pdb import Pdb, pdb_from_file
from kripo.protonate import protonate_protein, protonate_molecule, protonate_pdb


def test_protonate_protein(orig_pdb_3heg: Pdb):
    block = orig_pdb_3heg.model().chain('A').to_file_string('pdb')
    # Verify input contains no hydrogens
    assert len(orig_pdb_3heg.model().chain('A').atoms(element='H')) == 0

    hblock, err = protonate_protein(block)

    hpdb = pdb_string2pdb(hblock)
    assert len(hpdb.model().chain('A').atoms(element='H')) == 2702


def nr_hydrogens(mol: Mol):
    return mol.GetNumAtoms() - mol.GetNumHeavyAtoms()


def test_protonate_molecule(ligand_expo_dict_fixture):
    mol = ligand_expo_dict_fixture['3heg_BAX_1_A_1']
    # Verify input contains no hydrogens
    assert nr_hydrogens(mol) == 0, 'Contains no hydrogens before protonation'

    hmol = protonate_molecule(mol)

    assert nr_hydrogens(hmol) == 16


def test_protonate_pdb_3heg(orig_pdb_3heg: Pdb):
    protonated = protonate_pdb(orig_pdb_3heg)
    assert len(protonated.model().atoms(element='H')) == 2702


def test_protonate_pdb_3rze():
    unprotonated = pdb_from_file('tests/fixtures/3RZE.pdb', hydrogenate=False, clean=False)
    protonated = protonate_pdb(unprotonated)
    assert len(protonated.model().atoms(element='H')) == 3574


def test_protonate_pdb_5is0():
    unprotonated = pdb_from_file('tests/fixtures/5IS0.pdb', hydrogenate=False, clean=False)
    protonated = protonate_pdb(unprotonated)
    assert len(protonated.model().atoms(element='H')) == 11072


def test_protonate_pdb_2muv():
    unprotonated = pdb_from_file('tests/fixtures/2MUV.pdb', hydrogenate=False, clean=False)
    protonated = protonate_pdb(unprotonated)
    assert len(protonated.model().atoms(element='H')) == 982


def pdb_string2pdb(block):
    return pdb_dict_to_pdb(pdb_string_to_pdb_dict(block))
