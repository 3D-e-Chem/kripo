from kripo.fragment import Fragment
from kripo.pharmacophore import from_fragment, Feature


def test_from_fragment__fragment2_3heg_bax(fragment2_3heg_bax: Fragment):
    pharmacophore = from_fragment(fragment2_3heg_bax)

    expected_features = {
        Feature('NEGC', [0.5855, 4.674, 27.41]),
        Feature('LIPO', [-4.515, -1.871, 26.27]),
        Feature('LIPO', [-5.434, 0.1428, 24.48]),
        Feature('LIPO', [-3.722, -1.799, 26.57]),
        Feature('LIPO', [-1.103, 1.12, 22.3]),
        Feature('LIPO', [6.189, 3.114, 22.44]),
        Feature('HDON', [3.708, 4.227, 24.48]),
        Feature('HACC', [-1.68, 2.363, 22.7]),
        Feature('LIPO', [1.311, 0.2163, 23.44]),
        Feature('POSC', [-7.201, 1.861, 25.49]),
        Feature('LIPO', [-0.6538, 1.132, 21.84]),
        Feature('HACC', [-7.211, 1.939, 25.43]),
        Feature('LIPO', [-2.587, -1.334, 25.06]),
        Feature('LIPO', [-4.173, 3.063, 30.17]),
        Feature('LIPO', [3.625, 2.58, 25.24]),
        Feature('NEGC', [3.888, 3.075, 24.47]),
        Feature('LIPO', [-0.1571, 0.2585, 23.97]),
        Feature('NEGC', [2.105, 4.733, 25.37]),
    }

    assert pharmacophore.features == expected_features


def test_from_fragment__fragment1_3heg_bax(fragment1_3heg_bax):
    pharmacophore = from_fragment(fragment1_3heg_bax)

    for f in pharmacophore.features:
        p = f.position
        s = '{0} {1:.4} {2:.4} {3:.4}  0 0 0 0 0'.format(f.kind, p[0], p[1], p[2])
        print(s)

    expected = {
        Feature('LIPO', [6.069, 0.3048, 16.1]),
        Feature('LIPO', [-5.434, 0.1428, 24.48]),
        Feature('LIPO', [3.625, 2.58, 25.24]),
        Feature('LIPO', [6.189, 3.114, 22.44]),
        Feature('HACC', [-1.68, 2.363, 22.7]),
        Feature('LIPO', [3.687, -1.433, 20.74]),
        Feature('LIPO', [1.311, 0.2163, 23.44]),
        Feature('LIPO', [-0.1571, 0.2585, 23.97]),
        Feature('LIPO', [-2.587, -1.334, 25.06]),
        Feature('LIPO', [-4.515, -1.871, 26.27]),
        Feature('LIPO', [-3.722, -1.799, 26.57]),
        Feature('LIPO', [-4.173, 3.063, 30.17]),
        Feature('LIPO', [3.063, 1.939, 17.2]),
        Feature('POSC', [-7.201, 1.861, 25.49]),
        Feature('HACC', [-7.211, 1.939, 25.43]),
        Feature('LIPO', [-0.6538, 1.132, 21.84]),
        Feature('LIPO', [-1.103, 1.12, 22.3]),
        Feature('HACC', [3.318, -3.779, 16.54]),
        Feature('HDON', [2.973, -3.443, 13.16]),
        Feature('LIPO', [5.75, 0.6242, 12.85]),
        Feature('LIPO', [6.235, 0.8313, 13.51]),
        Feature('HDON', [3.708, 4.227, 24.48]),
        Feature('NEGC', [3.888, 3.075, 24.47]),
        Feature('NEGC', [0.5855, 4.674, 27.41]),
        Feature('NEGC', [2.105, 4.733, 25.37]),
    }
    assert pharmacophore.features == expected
