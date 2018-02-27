from rdkit.Chem.rdmolfiles import MolToSmiles

from kripo.generator import add_fragment2db, build_frag_id
from kripodb.db import FragmentsDb


def test_add_fragment2db(pdb_3heg, ligand_3heg_bax, fragment25_3heg_bax):
    pdb = pdb_3heg
    ligand = ligand_3heg_bax
    frag_nr = 25
    fragment = fragment25_3heg_bax
    fragment.name = '3heg_BAX_frag25'
    db = FragmentsDb(':memory:')

    add_fragment2db(pdb, ligand, frag_nr, fragment, db)

    row = db['3heg_BAX_frag25']
    assert MolToSmiles(row['mol']) == '[*]C(=O)NC'
    del row['mol']  # ignore molecule, assert based on smile above
    expected_row = {
        'atom_codes': 'N30,O32,C31',
        'ec_number': None,
        'frag_id': '3heg_BAX_frag25',
        'frag_nr': 25,
        'hash_code': 'c63c1cbc858154f00089fe78fdd4b2ba',
        'het_chain': 'A',
        'het_code': 'BAX',
        'het_seq_nr': 1,
        'nr_r_groups': 1,
        'pdb_code': '3heg',
        'pdb_title': None,
        'prot_chain': 'A',
        'prot_name': None,
        'rowid': 1,
        'smiles': '[*]C(=O)NC',
        'uniprot_acc': None,
        'uniprot_name': None,
    }
    assert row == expected_row


def test_build_frag_id(pdb_3heg, ligand_3heg_bax):

    frag_id = build_frag_id(pdb_3heg, ligand_3heg_bax, 13)

    assert frag_id == '3heg_BAX_frag13'
