from kripo.fragment import Fragment
from kripo.ligand import Ligand


def test_parent(ligand_3heg_bax: Ligand, fragment1_3heg_bax: Fragment):
    assert fragment1_3heg_bax.parent == ligand_3heg_bax.molecule


def test_atom_names__when_fragment_is_ligand(ligand_3heg_bax: Ligand, fragment1_3heg_bax: Fragment):
    names = fragment1_3heg_bax.atom_names()

    expected_names = [a.name() for a in ligand_3heg_bax.molecule.atoms()]
    assert len(names) > 0
    assert set(names) == set(expected_names)


def test_atom_names__when_fragment_is_not_ligand(ligand_3heg_bax: Ligand, fragment2_3heg_bax: Fragment):
    names = fragment2_3heg_bax.atom_names()

    expected_names = {a.name() for a in ligand_3heg_bax.molecule.atoms()}
    assert 0 < len(names) < len(expected_names)
    assert expected_names.issuperset(names)


def test_site(fragment2_3heg_bax: Fragment, ligand_3heg_bax: Ligand):
    site = fragment2_3heg_bax.site()

    assert site.ligand() == ligand_3heg_bax.molecule
    assert len(site.residues()) == 19


def test_site1(fragment1_3heg_bax: Fragment):
    site = fragment1_3heg_bax.site()
    seq_nrs = [int(r.residue_id().replace('A', '')) for r in site.residues()]
    seq_nrs.sort()
    [print(s) for s in seq_nrs]
    print(len(seq_nrs))
    assert False


def test_nr_r_groups(fragment2_3heg_bax: Fragment):
    nr = fragment2_3heg_bax.nr_r_groups()

    assert nr == 1


def test_hash_code__when_fragment_is_ligand(fragment1_3heg_bax: Fragment):
    hash_code = fragment1_3heg_bax.hash_code()

    expected_hash_code = '75add68f789c24266bce8e76d474acee'
    assert hash_code == expected_hash_code


def test_hash_code__when_fragment_is_not_ligand(fragment2_3heg_bax: Fragment):
    hash_code = fragment2_3heg_bax.hash_code()

    expected_hash_code = '930391721063a9b18730c8adfffef29a'
    assert hash_code == expected_hash_code


def test_mol_block__when_fragment_is_ligand(fragment1_3heg_bax: Fragment):
    mol_block = fragment1_3heg_bax.mol_block('3HEG_BAX_frag1')

    assert '3HEG_BAX_frag1' in mol_block
    assert 'V2000' in mol_block
    assert 'END' in mol_block
    assert len(mol_block.split('\n')) == 72
    # TODO compare whole block instead of pieces


def test_mol_block__when_fragment_is_not_ligand(fragment2_3heg_bax: Fragment):
    mol_block = fragment2_3heg_bax.mol_block('3HEG_BAX_frag2')

    assert '3HEG_BAX_frag2' in mol_block
    assert 'V2000' in mol_block
    assert 'END' in mol_block
    assert len(mol_block.split('\n')) == 51
    # TODO compare whole block instead of pieces
