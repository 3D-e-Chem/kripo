from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.structures import Model

from kripo.pdb import Pdb
from kripo.protonate import protonate_protein, protonate_ligand, protonate


def test_protonate_protein(pdb_3heg: Pdb):
        block = pdb_3heg.model.chain('A').to_file_string('pdb')

        hblock, err = protonate_protein(block)

        hpdb = pdb_string2pdb(hblock)
        assert len(hpdb.model().chain('A').atoms(element='H')) == 2724


def test_protonate_ligand(pdb_3heg: Pdb):
        # Create model which only contains ligand molecules
        ligands_model = Model()
        [ligands_model.add_molecule(l.molecule) for l in pdb_3heg.ligands()]
        block = ligands_model.to_file_string('pdb')
        # Verify input contains no hydrogens
        assert len(ligands_model.atoms(element='H')) == 0

        hblock = protonate_ligand(block)

        hpdb = pdb_string2pdb(hblock)
        hmodel = hpdb.model()
        assert len(hmodel.atoms(element='H')) == 16


def test_protonate_pdb_3heg(pdb_3heg: Pdb):
    protonated_model = protonate(pdb_3heg.model)
    assert len(protonated_model.atoms(element='H')) == 2740


def test_protonate_pdb_3rze(pdb_3rze: Pdb):
    protonated_model = protonate(pdb_3rze.model)
    assert len(protonated_model.atoms(element='H')) == 3642


def test_protonate_pdb_5is0(pdb_5is0: Pdb):
    protonated_model = protonate(pdb_5is0.model)
    assert len(protonated_model.atoms(element='H')) == 11012


def test_protonate_pdb_2muv(pdb_2muv: Pdb):
    protonated_model = protonate(pdb_2muv.model)
    protonated_model.save('2muv.protonated.pdb')
    assert len(protonated_model.atoms(element='H')) == 1069


def pdb_string2pdb(block):
    return pdb_dict_to_pdb(pdb_string_to_pdb_dict(block))
