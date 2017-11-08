from kripo.ligand import Ligand


def test_name(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.name() == 'BAX'


def test_pdb_block(ligand_3heg_bax: Ligand):
    block = ligand_3heg_bax.pdb_block()
    assert 'HETATM' in block
    assert 'BAX' in block
    assert 'CL11' in block
    assert 'CONECT' in block


def test_fragments(ligand_3heg_bax: Ligand):
    fragments = ligand_3heg_bax.fragments()

    assert len(fragments) == 10


def test_id(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.id() == 'A1'


def test_chain(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.chain() == 'A'


def test_seq_nr(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.seq_nr() == 1
