from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.structures import Model

from kripo.pdb import Pdb, pdb_from_file, ligands
from kripo.protonate import protonate_protein, protonate_ligand, protonate


def test_protonate_protein(pdb_3heg: Pdb):
        block = pdb_3heg.model().chain('A').to_file_string('pdb')

        hblock, err = protonate_protein(block)

        hpdb = pdb_string2pdb(hblock)
        assert len(hpdb.model().chain('A').atoms(element='H')) == 2724


def test_protonate_ligand():
    # Create model which only contains ligand molecules
    filename = 'tests/fixtures/3HEG.pdb'
    pdb_3heg = pdb_from_file(filename, hydrogenate=False, clean=False)
    ligands_model = Model()
    [ligands_model.add_molecule(l.molecule) for l in ligands(pdb_3heg)]
    block = ligands_model.to_file_string('pdb')
    # Verify input contains no hydrogens
    assert len(ligands_model.atoms(element='H')) == 0

    hblock = protonate_ligand(block)

    hpdb = pdb_string2pdb(hblock)
    hmodel = hpdb.model()
    assert len(hmodel.atoms(element='H')) == 18


def test_protonate_pdb_3heg():
    unprotonated = pdb_from_file('tests/fixtures/3HEG.pdb', hydrogenate=False, clean=False)
    protonated = protonate(unprotonated)
    assert len(protonated.model().atoms(element='H')) == 2958


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
