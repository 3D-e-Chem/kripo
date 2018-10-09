import pytest

from atomium.structures import Residue

from kripo.fragment import Fragment, is_residue_nearby, bounding_box, bounding_boxes_overlap
from kripo.ligand import Ligand
from rdkit.Chem.rdmolfiles import MolToSmiles


def test_parent(ligand_3heg_bax: Ligand, fragment1_3heg_bax: Fragment):
    assert fragment1_3heg_bax.parent == ligand_3heg_bax.atomium_mol


def test_atom_names__when_fragment_is_whole_ligand(fragment1_3heg_bax: Fragment):
    names = fragment1_3heg_bax.atom_names()

    expected = {
        'C1',
        'C13',
        'C16',
        'C17',
        'C18',
        'C19',
        'C2',
        'C20',
        'C21',
        'C23',
        'C24',
        'C25',
        'C27',
        'C28',
        'C29',
        'C3',
        'C31',
        'C4',
        'C5',
        'C6',
        'C7',
        'CL11',
        'F10',
        'F8',
        'F9',
        'N12',
        'N14',
        'N26',
        'N30',
        'O15',
        'O22',
        'O32',
    }

    assert set(names) == expected


def test_atom_names__when_fragment_is_part_ofligand(fragment2_3heg_bax: Fragment):
    names = set(fragment2_3heg_bax.atom_names())

    expected = {
        'C1',
        'C13',
        'C16',
        'C17',
        'C18',
        'C19',
        'C2',
        'C20',
        'C21',
        'C23',
        'C24',
        'C25',
        'C27',
        'C28',
        'C3',
        'C4',
        'C5',
        'C6',
        'C7',
        'CL11',
        'F10',
        'F8',
        'F9',
        'N12',
        'N14',
        'N26',
        'O15',
        'O22'
    }

    assert names == expected


def test_atom_names__of_small_fragment(fragment25_3heg_bax: Fragment):
    names = set(fragment25_3heg_bax.atom_names())

    expected = {'O32', 'N30', 'C29', 'C31'}

    assert names == expected


def seq_nrs_of_site(site, chain='A'):
    return {int(r.residue_id().replace(chain, '')) for r in site.residues()}


def test_site__fragment2_3heg_bax(fragment2_3heg_bax: Fragment):
    site = fragment2_3heg_bax.site()

    site_ligand = site.ligand()
    assert len(site_ligand.atoms()) == 28
    assert site_ligand.name() == 'fragment of BAX'
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

    expected_hash_code = '51980afadf63a7e7c8e95ba4f5e412ee'
    assert hash_code == expected_hash_code


def test_hash_code__when_fragment_is_not_ligand(fragment2_3heg_bax: Fragment):
    hash_code = fragment2_3heg_bax.hash_code()

    expected_hash_code = '88d5b594a6326908206691a497bdbe19'
    assert hash_code == expected_hash_code


def test_set_name(fragment1_3heg_bax: Fragment):
    name = '3HEG_BAX_frag1'
    fragment1_3heg_bax.name = name

    assert fragment1_3heg_bax.name == name
    assert fragment1_3heg_bax.molecule.GetProp('_Name') == name
    fragment1_3heg_bax.name = ''


def test_get_name(fragment1_3heg_bax: Fragment):
    assert fragment1_3heg_bax.name == ''


def test_is_residue_nearby__nothing_in__false():
    fragment_atoms = []
    residue = Residue()
    radius = 6.0
    with pytest.raises(ValueError):
        is_residue_nearby(fragment_atoms, residue, radius)


def test_smiles(fragment1_3heg_bax: Fragment):
    mol = fragment1_3heg_bax.smiles()

    expected = '[H]c1nc(C(=O)N([H])C([H])([H])[H])c([H])c(Oc2c([H])c([H])c(N([H])C(=O)N([H])c3c([H])c([H])c(Cl)c(C(F)' \
               '(F)F)c3[H])c([H])c2[H])c1[H]'
    assert mol == expected


def test_unprotonated_molecule(fragment1_3heg_bax: Fragment):
    mol = fragment1_3heg_bax.unprotonated_molecule()

    expected = 'CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1'
    assert MolToSmiles(mol) == expected


@pytest.mark.parametrize(('res', 'expected_overlaps'), (
    (((-1, 0), (-1, 0), (-1, 0)), True),  # same
    (((-.75, -0.25), (-.75, -.25), (-.75, -0.25)), True),  # inside
    (((-3, 2), (-3, 2), (-3, 2)), True),  # inside radius
    (((-10, 9), (-10, 9), (-10, 9)), True),  # frag inside res
    (((6, 7), (-1, 0), (-1, 0)), False),  # side x
    (((-1, 0), (6, 7), (-1, 0)), False),  # side y
    (((-1, 0), (-1, 0), (6, 7)), False),  # side z
    (((-6, -7), (-1, 0), (-1, 0)), False),  # side -x
    (((-1, 0), (-6, -7), (-1, 0)), False),  # side -y
    (((-1, 0), (-1, 0), (-6, -7)), False),  # side -z
    (((3, 6), (3, 6), (3, 6)), True),  # corner
    (((3, 6), (3, 6), (-7, -4)), True),  # corner
    (((3, 6), (-7, -4), (3, 6)), True),  # corner
    (((-7, -4), (3, 6), (3, 6)), True),  # corner
    (((-7, -4), (-7, -4), (3, 6)), True),  # corner
    (((-7, -4), (-7, -4), (-7, -4)), True),  # corner
    (((3, 6), (-7, -4), (-7, -4)), True),  # corner
    (((-7, -4), (3, 6), (-7, -4)), True),  # corner
    (((-2, 1), (3, 6), (3, 6)), True),  # edge
    (((3, 6), (-2, 1), (3, 6)), True),  # edge
    (((3, 6), (3, 6), (-2, 1)), True),  # edge
    (((-2, 1), (-7, -4), (-7, -4)), True),  # edge
    (((-7, -4), (-2, 1), (-7, -4)), True),  # edge
    (((-7, -4), (-7, -4), (-2, 1)), True),  # edge
    (((-2, 1), (-2, 1), (-7, -4)), True),  # side
    (((-2, 1), (-7, -4), (-2, 1)), True),  # side
    (((-7, -4), (-2, 1), (-2, 1)), True),  # side
    (((-17, -14), (-2, 1), (-2, 1)), False),  # out side
))
def test_bounding_boxes_overlap(res, expected_overlaps):
    frag = ((-1, 0), (-1, 0), (-1, 0))

    overlaps = bounding_boxes_overlap(frag, 5.0, res)
    assert overlaps == expected_overlaps


@pytest.mark.parametrize(('atoms', 'expected_bb',), (
    (
        (
            (0.0, 0.0, 0.0),
        ),
        (
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ),
    ), (
        (
            (0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0),
        ),
        (
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
        ),
    ), (
        (
            (1.0, 1.0, 1.0),
            (0.0, 0.0, 0.0),
        ),
        (
            (0.0, 1.0),
            (0.0, 1.0),
            (0.0, 1.0),
        ),
    ),
))
def test_bounding_box(atoms, expected_bb):
    bb = bounding_box(atoms)

    assert bb == expected_bb


def test_bounding_box_zero_atoms_value_error():
    atoms = list()

    with pytest.raises(ValueError):
        bounding_box(atoms)
