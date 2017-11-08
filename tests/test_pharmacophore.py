from kripo.fragment import Fragment
from kripo.pharmacophore import from_fragment, Feature


def test_from_fragment(fragment2_3heg_bax: Fragment):
    pharmacophore = from_fragment(fragment2_3heg_bax)

    expected_features = {
        Feature('LIPO', [-1.216, 1.885, 28.69]),
        Feature('LIPO', [-3.466, -2.357, 27.65]),
        Feature('AROM', [-0.8878, 1.379, 19.71]),
        Feature('HDON', [-3.334, -2.499, 24.26]),
        Feature('LIPO', [-7.188, 0.2354, 24.1]),
        Feature('HDON', [-0.03166, 4.629, 22.5]),
        Feature('LIPO', [1.948, -0.9389, 22.22]),
        Feature('LIPO', [-0.9244, -0.6887, 25.29]),
        Feature('LIPO', [-7.558, -1.195, 27.19]),
        Feature('NEGC', [3.888, 3.075, 24.47]),
        Feature('HDON', [3.708, 4.227, 24.48]),
        Feature('HDON', [1.149, 5.466, 26.76]),
        Feature('NEGC', [0.5654, 6.421, 26.43]),
        Feature('NEGC', [3.868, 4.822, 23.49]),
        Feature('NEGC', [0.5855, 4.674, 27.41]),
        Feature('NEGC', [2.105, 4.733, 25.37]),
        Feature('LIPO', [-5.957, 4.507, 30.93]),
    }
    assert pharmacophore.features == expected_features


def test_from_fragment1(fragment1_3heg_bax):
    pharmacophore = from_fragment(fragment1_3heg_bax)

    for f in pharmacophore.features:
        p = f.position
        s = '{0} {1:.4} {2:.4} {3:.4}  0 0 0 0 0'.format(f.kind, p[0], p[1], p[2])
        print(s)

    expected = {
        Feature('HACC', [-1.7076, 2.2682, 22.7126]),
        Feature('HDON', [-0.0317, 4.6294, 22.4973]),
        Feature('HDON', [3.9657, -4.2182, 19.4535]),
        Feature('AROM', [-4.7420, 5.8751, 25.9774]),
        Feature('HDON', [3.7079, 4.2267, 24.4837]),
        Feature('NEGC', [3.8882, 3.0747, 24.4667]),
        Feature('HDON', [1.1487, 5.4662, 26.7621]),
        Feature('NEGC', [0.5654, 6.4209, 26.4331]),
        Feature('NEGC', [0.5855, 4.6737, 27.4061]),
        Feature('NEGC', [2.1046, 4.7328, 25.3711]),
        Feature('LIPO', [-0.9244, -0.6887, 25.2947]),
        Feature('LIPO', [1.3499, 0.1450, 23.3679]),
        Feature('LIPO', [-0.2041, 0.2003, 24.0542]),
        Feature('LIPO', [-2.4841, -1.2942, 25.0740]),
        Feature('LIPO', [-5.5423, 0.1490, 24.4589]),
        Feature('LIPO', [-2.2450, 5.3313, 29.7474]),
        Feature('LIPO', [-2.8484, 1.4969, 28.8643]),
        Feature('LIPO', [-0.7917, 3.4229, 28.1349]),
        Feature('LIPO', [3.7043, 2.5349, 25.3078]),
        Feature('LIPO', [1.4568, 1.3568, 25.5041]),
        Feature('LIPO', [5.7106, 2.4540, 25.7827]),
        Feature('LIPO', [5.4877, -3.2515, 15.6444]),
        Feature('LIPO', [-0.1181, 1.2310, 19.4872]),
        Feature('LIPO', [-0.7460, 1.0889, 21.8024]),
        Feature('LIPO', [6.2402, 3.1961, 22.3875]),
        Feature('LIPO', [-4.3778, 0.0278, 30.1472]),
        Feature('HACC', [3.3976, -3.8262, 16.5781]),
        Feature('HDON', [2.9732, -3.4432, 13.1614]),
        Feature('LIPO', [1.9186, -2.1014, 15.8518]),
        Feature('LIPO', [1.0415, -0.2159, 14.6303]),
        Feature('LIPO', [2.9680, 1.9694, 17.1866]),
        Feature('LIPO', [4.0759, -3.1885, 20.6854]),
        Feature('LIPO', [3.7108, -1.5415, 20.7330]),
        Feature('LIPO', [5.8368, 0.5870, 12.9128]),
        Feature('LIPO', [5.5328, 2.5229, 16.8755]),
        Feature('LIPO', [6.1775, 0.2833, 16.1111]),
        Feature('LIPO', [-4.4496, -1.8999, 26.3564]),
    }
    assert pharmacophore.features == expected
