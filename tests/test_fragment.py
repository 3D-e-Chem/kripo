import pytest

from atomium.structures import Residue

from kripo.fragment import Fragment, is_residue_nearby
from kripo.ligand import Ligand


def test_parent(ligand_3heg_bax: Ligand, fragment1_3heg_bax: Fragment):
    assert fragment1_3heg_bax.parent == ligand_3heg_bax.molecule


def test_atom_names__when_fragment_is_ligand(ligand_3heg_bax: Ligand, fragment1_3heg_bax: Fragment):
    names = fragment1_3heg_bax.atom_names()

    expected_names = {a.name() for a in ligand_3heg_bax.molecule.atoms() if a.name() != 'H'}
    assert len(names) > 0
    assert set(names) == expected_names


def test_atom_names__when_fragment_is_not_ligand(ligand_3heg_bax: Ligand, fragment2_3heg_bax: Fragment):
    names = fragment2_3heg_bax.atom_names()

    expected_names = {a.name() for a in ligand_3heg_bax.molecule.atoms()}
    expected_names.remove('H')
    assert 0 < len(names) < len(expected_names)
    assert expected_names.issuperset(names)


def test_site__fragment2_3heg_bax(fragment2_3heg_bax: Fragment):
    site = fragment2_3heg_bax.site()

    site_ligand = site.ligand()
    assert len(site_ligand.atoms()) == 20
    assert site_ligand.name() == 'BAX'
    assert len(site.residues()) == 25


def test_site__fragment1_3heg_bax(fragment1_3heg_bax: Fragment):
    site = fragment1_3heg_bax.site()

    seq_nrs = [int(r.residue_id().replace('A', '')) for r in site.residues()]
    seq_nrs.sort()
    expected = {138, 140, 141, 146, 147, 148, 149, 151, 157, 30, 35, 166, 167, 40, 169, 168, 38, 51, 53, 71, 74, 75, 78, 83, 84, 86, 104, 106, 107, 108, 109, 110}
    assert expected == set(seq_nrs)


def test_site_yasara_fragment(yasara_fragment2_3heg_bax: Fragment):
    site = yasara_fragment2_3heg_bax.site()

    seq_nrs = {int(r.residue_id().replace('A', '')) for r in site.residues()}
    filename = 'tests/fixtures/3HEG.frag2.site.pdb'
    from atomium.files import pdb_from_file
    expected_pdb = pdb_from_file(filename).model()
    expected = {int(r.residue_id().replace('A', '')) for r in expected_pdb.residues()}
    # The yasara selection and this selection implementation have different prep steps and behaviors
    # They should show enough overlap
    assert len(seq_nrs & expected) > 20


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


def test_is_residue_nearby__nothing_in__false():
    fragment_atoms = set()
    residue = Residue()
    radius = 6.0
    assert not is_residue_nearby(fragment_atoms, residue, radius)


def test_is_residue_nearby__radiustoobig_valueerror():
    fragment_atoms = set()
    residue = Residue()
    radius = 999999999.0
    with pytest.raises(ValueError) as excinfo:
        is_residue_nearby(fragment_atoms, residue, radius)

    assert 'Radius must be smaller than' in str(excinfo.value)
