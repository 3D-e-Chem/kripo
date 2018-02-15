import pytest

from kripo.fingerprint.utils import fuzzy_offsets


@pytest.mark.parametrize('factor,shape,expected', (
    (0, 'all', set()),
    (0, 'one', set()),
    (0, 'v1', set()),
    (1, 'all', set()),
    (1, 'one', set()),
    (1, 'v1', set()),
))
def test_fuzzy_offsets(factor, shape, expected):
    result = set(fuzzy_offsets(factor, shape))

    assert result == expected
