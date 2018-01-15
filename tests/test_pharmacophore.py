from kripo.fragment import Fragment
from kripo.pharmacophore import from_fragment, Feature


def test_from_fragment__fragment2_3heg_bax(fragment2_3heg_bax: Fragment):
    pharmacophore = from_fragment(fragment2_3heg_bax)

    expected_features = {
        Feature('LIPO', [6.069, 0.3048, 16.1]),
        Feature('LIPO', [5.448, 2.463, 16.84]),
        Feature('HACC', [3.966, -4.218, 19.45]),
        Feature('AROM', [-4.742, 5.875, 25.98]),
        Feature('POSC', [-7.201, 1.861, 25.49]),
        Feature('HDON', [-7.211, 1.939, 25.43]),
        Feature('LIPO', [-0.01681, 1.203, 19.52]),
        Feature('LIPO', [-0.6538, 1.132, 21.84]),
        Feature('HDON', [3.318, -3.779, 16.54]),
        Feature('LIPO', [1.96, -2.005, 15.89]),
        Feature('LIPO', [-2.587, -1.334, 25.06]),
        Feature('LIPO', [-0.9244, -0.6887, 25.29]),
        Feature('LIPO', [1.311, 0.2163, 23.44]),
        Feature('LIPO', [-0.1571, 0.2585, 23.97]),
        Feature('LIPO', [6.189, 3.114, 22.44]),
        Feature('LIPO', [-4.964, -0.1787, 29.82]),
        Feature('HACC', [3.708, 4.227, 24.48]),
        Feature('HACC', [1.149, 5.466, 26.76]),
        Feature('NEGC', [3.888, 3.075, 24.47]),
        Feature('NEGC', [0.5855, 4.674, 27.41]),
        Feature('NEGC', [2.105, 4.733, 25.37]),
        Feature('NEGC', [0.5654, 6.421, 26.43]),
        Feature('LIPO', [3.625, 2.58, 25.24]),
        Feature('LIPO', [1.454, 1.444, 25.43]),
        Feature('LIPO', [-5.434, 0.1428, 24.48]),
        Feature('HDON', [-1.68, 2.363, 22.7]),
        Feature('LIPO', [4.076, -3.188, 20.69]),
        Feature('LIPO', [3.687, -1.433, 20.74]),
        Feature('HACC', [-0.03166, 4.629, 22.5]),
        Feature('LIPO', [-4.515, -1.871, 26.27]),
        Feature('LIPO', [-4.173, 3.063, 30.17]),
        Feature('LIPO', [-2.182, 5.39, 29.68]),
        Feature('LIPO', [-2.955, 1.471, 28.88]),
        Feature('LIPO', [-0.7636, 3.524, 28.1]),
        Feature('LIPO', [3.063, 1.939, 17.2]),
        Feature('LIPO', [5.668, 2.359, 25.82]),
        Feature('LIPO', [1.213, -0.4437, 14.56]),
        Feature('LIPO', [5.396, -3.238, 15.71]),
        Feature('HACC', [2.973, -3.443, 13.16]),
        Feature('LIPO', [5.585, -0.8829, 10.59])
    }

    assert pharmacophore.features == expected_features


def test_from_fragment__fragment1_3heg_bax(fragment1_3heg_bax):
    pharmacophore = from_fragment(fragment1_3heg_bax)

    for f in pharmacophore.features:
        p = f.position
        s = '{0} {1:.4} {2:.4} {3:.4}  0 0 0 0 0'.format(f.kind, p[0], p[1], p[2])
        print(s)

    expected = {
        Feature('LIPO', [3.625, 2.58, 25.24]),
        Feature('LIPO', [1.454, 1.444, 25.43]),
        Feature('LIPO', [-5.434, 0.1428, 24.48]),
        Feature('HACC', [3.708, 4.227, 24.48]),
        Feature('NEGC', [3.888, 3.075, 24.47]),
        Feature('HACC', [1.149, 5.466, 26.76]),
        Feature('NEGC', [0.5855, 4.674, 27.41]),
        Feature('NEGC', [2.105, 4.733, 25.37]),
        Feature('NEGC', [0.5654, 6.421, 26.43]),
        Feature('LIPO', [5.396, -3.238, 15.71]),
        Feature('LIPO', [6.189, 3.114, 22.44]),
        Feature('LIPO', [5.668, 2.359, 25.82]),
        Feature('LIPO', [4.076, -3.188, 20.69]),
        Feature('HDON', [-1.68, 2.363, 22.7]),
        Feature('LIPO', [3.687, -1.433, 20.74]),
        Feature('HACC', [-0.03166, 4.629, 22.5]),
        Feature('LIPO', [-4.515, -1.871, 26.27]),
        Feature('LIPO', [5.75, 0.6242, 12.85]),
        Feature('LIPO', [5.585, -0.8829, 10.59]),
        Feature('LIPO', [-4.173, 3.063, 30.17]),
        Feature('LIPO', [-2.182, 5.39, 29.68]),
        Feature('LIPO', [-0.7636, 3.524, 28.1]),
        Feature('LIPO', [-2.955, 1.471, 28.88]),
        Feature('LIPO', [3.063, 1.939, 17.2]),
        Feature('LIPO', [6.069, 0.3048, 16.1]),
        Feature('LIPO', [5.448, 2.463, 16.84]),
        Feature('AROM', [-4.742, 5.875, 25.98]),
        Feature('POSC', [-7.201, 1.861, 25.49]),
        Feature('HDON', [-7.211, 1.939, 25.43]),
        Feature('HACC', [3.966, -4.218, 19.45]),
        Feature('LIPO', [-0.01681, 1.203, 19.52]),
        Feature('LIPO', [-0.6538, 1.132, 21.84]),
        Feature('HDON', [3.318, -3.779, 16.54]),
        Feature('HACC', [2.973, -3.443, 13.16]),
        Feature('LIPO', [1.96, -2.005, 15.89]),
        Feature('LIPO', [1.213, -0.4437, 14.56]),
        Feature('LIPO', [-4.964, -0.1787, 29.82]),
        Feature('LIPO', [-0.9244, -0.6887, 25.29]),
        Feature('LIPO', [1.311, 0.2163, 23.44]),
        Feature('LIPO', [-0.1571, 0.2585, 23.97]),
        Feature('LIPO', [-2.587, -1.334, 25.06]),
    }
    assert pharmacophore.features == expected
