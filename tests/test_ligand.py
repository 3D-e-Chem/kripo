from kripo.ligand import Ligand
from rdkit.Chem.rdmolfiles import  MolToSmiles


def test_name(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.name() == 'BAX'


def test_fragments(ligand_3heg_bax: Ligand):
    fragments = ligand_3heg_bax.fragments()

    assert len(fragments) == 28


def test_id(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.id() == 'A1'


def test_chain(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.chain() == 'A'


def test_seq_nr(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.seq_nr() == 1


def test_rdkit_mol(ligand_3heg_bax: Ligand):
    mol = ligand_3heg_bax.rdkit_mol

    smiles = MolToSmiles(mol)
    assert 'c' in smiles  # aromatic bonds
    assert 'H' in smiles  # has hydrogens
    assert '=' in smiles  # double bonds
