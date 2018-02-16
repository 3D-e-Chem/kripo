import pytest

from kripo.fingerprint.utils import fuzzy_offsets


@pytest.mark.parametrize('factor,shape,expected', (
    (0, 'all', set()),
    (0, 'one', set()),
    (0, 'v1', set()),
    (1, 'all', {
        (0, 0, -1),
        (0, 0, 1),
        (0, -1, 0),
        (0, 1, 0),
        (0, 1, -1),
        (0, 1, 1),
        (0, -1, -1),
        (0, -1, 1),
        (-1, 0, 0),
        (1, 0, 0),
        (-1, 0, 1),
        (-1, 0, -1),
        (-1, 1, 0),
        (-1, -1, 0),
        (-1, 1, -1),
        (-1, 1, 1),
        (-1, -1, -1),
        (-1, -1, 1),
        (1, 0, 1),
        (1, 0, -1),
        (1, 1, 0),
        (1, -1, 0),
        (1, 1, -1),
        (1, 1, 1),
        (1, -1, -1),
        (1, -1, 1),
    }),
    (1, 'one', {
        (0, 0, -1), (0, 0, 1),
        (0, -1, 0), (0, 1, 0),
        (-1, 0, 0), (1, 0, 0),
    }),
    (1, 'v1', {
        (1, 2, 1), (1, 1, -1), (1, -1, -1),
        (1, 0, -1), (0, -1, -1), (-1, -1, -1),
        (1, 2, 0)
    }),
    (2, 'one', {
        (0, 0, -1), (0, 0, 1), (0, -1, 0),
        (0, 1, 0), (-1, 0, 0), (1, 0, 0),
        (0, 0, -2), (0, 0, 2), (0, -2, 0),
        (0, 2, 0), (-2, 0, 0), (2, 0, 0),
    }),
))
def test_fuzzy_offsets(factor, shape, expected):
    result = set(fuzzy_offsets(factor, shape))

    assert result == expected


def test_fuzzy_offsets__wrongfactor_valueerror():
    with pytest.raises(ValueError):
        list(fuzzy_offsets(-1, 'all'))


def test_fuzzy_offsets__wrongshape_valueerror():
    with pytest.raises(ValueError):
        list(fuzzy_offsets(1, 'peer'))
