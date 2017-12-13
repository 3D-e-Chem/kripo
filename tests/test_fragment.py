import pytest

from atomium.structures import Residue

from kripo.fragment import Fragment, is_residue_nearby
from kripo.ligand import Ligand


def test_parent(ligand_3heg_bax: Ligand, fragment1_3heg_bax: Fragment):
    assert fragment1_3heg_bax.parent == ligand_3heg_bax.molecule


def test_atom_names__when_fragment_is_whole_ligand(ligand_3heg_bax: Ligand, fragment1_3heg_bax: Fragment):
    names = fragment1_3heg_bax.atom_names()

    expected_names = {a.name() for a in ligand_3heg_bax.molecule.atoms()}
    expected_names.remove('H')
    assert len(names) > 0
    assert set(names) == expected_names


def test_atom_names__when_fragment_is_part_ofligand(ligand_3heg_bax: Ligand, fragment2_3heg_bax: Fragment):
    names = fragment2_3heg_bax.atom_names()

    expected_names = {a.name() for a in ligand_3heg_bax.molecule.atoms()}
    expected_names.remove('H')
    assert 0 < len(names) < len(expected_names)
    assert expected_names.issuperset(names)


def seq_nrs_of_site(site):
    return {int(r.residue_id().replace('A', '')) for r in site.residues()}


def test_site__fragment2_3heg_bax(fragment2_3heg_bax: Fragment):
    site = fragment2_3heg_bax.site()

    site_ligand = site.ligand()
    assert len(site_ligand.atoms()) == 39
    assert site_ligand.name() == 'BAX'
    seq_nrs = seq_nrs_of_site(site)
    expected = {
        30,
        35,
        36,
        38,
        51,
        53,
        55,
        70,
        71,
        74,
        75,
        78,
        83,
        84,
        86,
        104,
        106,
        107,
        108,
        109,
        110,
        138,
        140,
        141,
        146,
        147,
        148,
        149,
        151,
        155,
        157,
        166,
        167,
        168,
        169}
    assert expected == seq_nrs


def test_site__fragment1_3heg_bax(fragment1_3heg_bax: Fragment):
    site = fragment1_3heg_bax.site()

    seq_nrs = seq_nrs_of_site(site)

    """Expected seq nrs was calculated using Yasara script::

        LoadPDB tests/fixtures/3HEG.prepped.pdb,Center=No,Correct=No
        DelRes HOH
        DelRes protein with distance > 6 from bax
        DelRes BAX
        SavePDB 1,tests/fixtures/3heg_BAX_frag1.prepped.site.pdb,Format=PDB,Transform=No

    and the Python snippet::

        from atomium.files import pdb_from_file
        p = pdb_from_file('3heg_BAX_frag1.prepped.site.pdb').model()
        seq_nrs = sorted([int(r.residue_id().replace('A', '')) for r in p.residues()])
    """
    expected = {
        30,
        35,
        36,
        38,
        40,
        51,
        53,
        55,
        70,
        71,
        74,
        75,
        78,
        83,
        84,
        86,
        104,
        106,
        107,
        108,
        109,
        110,
        111,
        138,
        140,
        141,
        146,
        147,
        148,
        149,
        151,
        155,
        157,
        166,
        167,
        168,
        169
    }
    assert expected == seq_nrs


def test_nr_r_groups(fragment2_3heg_bax: Fragment):
    nr = fragment2_3heg_bax.nr_r_groups()

    assert nr == 1


def test_hash_code__when_fragment_is_ligand(fragment1_3heg_bax: Fragment):
    hash_code = fragment1_3heg_bax.hash_code()

    expected_hash_code = '75add68f789c24266bce8e76d474acee'
    assert hash_code == expected_hash_code


def test_hash_code__when_fragment_is_not_ligand(fragment2_3heg_bax: Fragment):
    hash_code = fragment2_3heg_bax.hash_code()

    expected_hash_code = 'cf6e9dedd21a2f6597e898664ec51977'
    assert hash_code == expected_hash_code


def test_set_name(fragment1_3heg_bax: Fragment):
    name = '3HEG_BAX_frag1'
    fragment1_3heg_bax.name = name

    assert fragment1_3heg_bax.name == name
    assert fragment1_3heg_bax.molecule.GetProp('_Name') == name


def test_get_name(fragment1_3heg_bax: Fragment):
    assert fragment1_3heg_bax.name == ''


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
