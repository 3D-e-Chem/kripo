from typing import Set
import json

from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.structures.chains import Site
import pytest

from kripo.ligand import Ligand
from kripo.pdb import pdb_from_file, ligands
from kripo.pharmacophore import Feature, from_site
from kripodb.pharmacophores import as_phar


def prep_site(block, ligand):
    model = pdb_dict_to_pdb(pdb_string_to_pdb_dict(block)).model()
    return Site(model, ligand=ligand)


def assert_features(expected: Set[Feature], result: Set[Feature]):
    sorted_expected = sorted(list(expected), key=lambda d: repr(d))
    sorted_result = sorted(list(result), key=lambda d: repr(d))
    assert sorted_expected == sorted_result


ALA = pytest.param(
    """ATOM    260  N   ALA A  40      10.884  -2.262  13.847  1.00 39.68           N  
ATOM    261  CA  ALA A  40      10.685  -3.687  14.015  1.00 40.02           C  
ATOM    262  C   ALA A  40      12.080  -4.263  14.391  1.00 39.68           C  
ATOM    263  O   ALA A  40      13.103  -3.778  13.914  1.00 38.58           O  
ATOM    264  CB  ALA A  40      10.148  -4.289  12.730  1.00 39.47           C  
ATOM         H   ALA A  40      10.831  -1.937  12.903  1.00 39.68           H  
ATOM         HA  ALA A  40       9.948  -3.922  14.797  1.00 40.02           H  
ATOM         HB1 ALA A  40      10.000  -5.371  12.865  1.00 39.47           H  
ATOM         HB2 ALA A  40       9.188  -3.816  12.476  1.00 39.47           H  
ATOM         HB3 ALA A  40      10.867  -4.117  11.915  1.00 39.47           H  
""",
    {
        Feature('HDON', [10.82, -1.872, 12.71]),
        Feature('HACC', [13.77, -3.462, 13.6]),
        Feature('LIPO', [11.39, -3.992, 11.32]),
        Feature('LIPO', [8.49, -3.472, 12.29]),
        Feature('LIPO', [9.865, -4.606, 12.05]),
        Feature('LIPO', [9.892, -6.158, 12.96]),
    },
    id='ALA'
)

ARG = pytest.param(
    """ATOM    508  N   ARG A  70       3.911   9.033  32.523  1.00 21.83           N  
ATOM    509  CA  ARG A  70       2.537   8.524  32.467  1.00 22.75           C  
ATOM    510  C   ARG A  70       2.491   7.281  31.559  1.00 22.34           C  
ATOM    511  O   ARG A  70       1.844   6.298  31.886  1.00 21.11           O  
ATOM    512  CB  ARG A  70       1.572   9.619  31.947  1.00 22.96           C  
ATOM    513  CG  ARG A  70       0.160   9.125  31.510  1.00 25.48           C  
ATOM    514  CD  ARG A  70      -0.680  10.221  30.789  1.00 23.42           C  
ATOM    515  NE  ARG A  70      -0.825  11.366  31.681  1.00 23.96           N  
ATOM    516  CZ  ARG A  70      -0.399  12.599  31.414  1.00 24.08           C  
ATOM    517  NH1 ARG A  70       0.143  12.900  30.234  1.00 20.76           N  
ATOM    518  NH2 ARG A  70      -0.560  13.537  32.331  1.00 24.15           N  
ATOM         HG2 ARG A  70       0.273   8.260  30.840  1.00 25.48           H  
ATOM         H   ARG A  70       4.022   9.955  32.153  1.00 21.83           H  
ATOM         HA  ARG A  70       2.212   8.242  33.480  1.00 22.75           H  
ATOM         HB2 ARG A  70       1.447  10.376  32.735  1.00 22.96           H  
ATOM         HB3 ARG A  70       2.045  10.122  31.091  1.00 22.96           H  
ATOM         HG3 ARG A  70      -0.389   8.776  32.397  1.00 25.48           H  
ATOM         HD3 ARG A  70      -1.669   9.824  30.515  1.00 23.42           H  
ATOM         HD2 ARG A  70      -0.186  10.526  29.855  1.00 23.42           H  
ATOM        HH11 ARG A  70       0.235  12.194  29.532  1.00 20.76           H  
ATOM         HE  ARG A  70      -1.278  11.214  32.559  1.00 23.96           H  
ATOM        HH21 ARG A  70      -0.994  13.314  33.204  1.00 24.15           H  
ATOM        HH12 ARG A  70       0.457  13.832  30.051  1.00 20.76           H  
ATOM        HH22 ARG A  70      -0.248  14.470  32.152  1.00 24.15           H  
""",
    {
        Feature('HDON', [-0.1856, 14.66, 32.12]), Feature('POSC', [0.265, 12.97, 29.97]),
        Feature('POSC', [-0.921, 11.09, 31.74]), Feature('LIPO', [1.356, 10.93, 33.31]),
        Feature('HDON', [4.044, 10.14, 32.08]), Feature('HACC', [1.42, 5.654, 32.1]),
        Feature('LIPO', [-0.7883, 8.522, 33.04]), Feature('HDON', [0.2534, 12.05, 29.39]),
        Feature('LIPO', [2.389, 10.49, 30.47]), Feature('LIPO', [0.3552, 7.631, 30.35]),
        Feature('HDON', [-1.369, 11.18, 32.73]), Feature('HDON', [-1.081, 13.27, 33.38]),
        Feature('POSC', [-0.5965, 13.75, 32.54]), Feature('HDON', [0.5198, 14.02, 30.01])
    },
    id='ARG'
)

ASN = pytest.param(
    """ATOM   1208  N   ASN A 155      -6.377   2.142  16.593  1.00 28.03           N  
ATOM   1209  CA  ASN A 155      -5.764   1.837  17.905  1.00 26.71           C  
ATOM   1210  C   ASN A 155      -5.967   0.396  18.367  1.00 26.62           C  
ATOM   1211  O   ASN A 155      -5.875   0.098  19.560  1.00 25.11           O  
ATOM   1212  CB  ASN A 155      -6.217   2.826  19.004  1.00 26.34           C  
ATOM   1213  CG  ASN A 155      -5.444   4.148  18.970  1.00 26.95           C  
ATOM   1214  OD1 ASN A 155      -6.013   5.220  19.149  1.00 28.68           O  
ATOM   1215  ND2 ASN A 155      -4.167   4.076  18.714  1.00 20.27           N  
ATOM         H   ASN A 155      -7.361   2.317  16.617  1.00 28.03           H  
ATOM         HA  ASN A 155      -4.683   1.963  17.742  1.00 26.71           H  
ATOM         HB2 ASN A 155      -7.291   3.032  18.885  1.00 26.34           H  
ATOM         HB3 ASN A 155      -6.087   2.356  19.990  1.00 26.34           H  
ATOM        HD21 ASN A 155      -3.735   3.185  18.573  1.00 20.27           H  
ATOM        HD22 ASN A 155      -3.620   4.911  18.659  1.00 20.27           H  
""",
    {
        Feature('HDON', [-3.649, 3.007, 18.54]), Feature('HACC', [-6.384, 5.919, 19.27]),
        Feature('HDON', [-3.511, 5.078, 18.65]), Feature('HDON', [-7.558, 2.352, 16.62]),
        Feature('HACC', [-5.815, -0.09534, 20.33])
    },
    id='ASN'
)

ASP = pytest.param(
    """ATOM   1309  N   ASP A 168      -2.779   2.584  22.274  1.00 25.68           N  
ATOM   1310  CA  ASP A 168      -3.014   3.950  21.910  1.00 28.49           C  
ATOM   1311  C   ASP A 168      -1.696   4.460  21.328  1.00 31.81           C  
ATOM   1312  O   ASP A 168      -0.684   4.563  22.039  1.00 32.09           O  
ATOM   1313  CB  ASP A 168      -3.390   4.777  23.160  1.00 28.61           C  
ATOM   1314  CG  ASP A 168      -3.835   6.192  22.812  1.00  28.1           C  
ATOM   1315  OD1 ASP A 168      -3.885   6.544  21.614  1.00 26.84           O  
ATOM   1316  OD2 ASP A 168      -4.158   6.949  23.725  1.00 25.16           O  
ATOM         H   ASP A 168      -1.863   2.400  22.631  1.00 25.68           H  
ATOM         HA  ASP A 168      -3.841   4.039  21.190  1.00 28.49           H  
ATOM         HB2 ASP A 168      -4.198   4.266  23.704  1.00 28.61           H  
ATOM         HB3 ASP A 168      -2.525   4.825  23.838  1.00 28.61           H  
""",
    {
        Feature('NEGC', [-2.999, 7.27, 20.33]), Feature('NEGC', [-3.584, 8.143, 24.82]),
        Feature('NEGC', [-4.244, 7.492, 22.49]), Feature('NEGC', [-5.467, 7.48, 24.71]),
        Feature('HACC', [-4.368, 7.442, 24.32]), Feature('NEGC', [-4.883, 6.607, 20.21]),
        Feature('HACC', [-0.03166, 4.629, 22.5]), Feature('HDON', [-1.68, 2.363, 22.7]),
        Feature('HACC', [-3.917, 6.769, 20.85])
    },
    id='ASP'
)

CYS = pytest.param(
    # The closest CYS to the ligand is 162, but it is more than 10 Angstroms away
    # Making difficult to view in the KNIME molviewer
    # The CYS 162 was translated 5, 5, 5 to bring it closer to the ligand
    # atomium_model.translate(5, 5, 5)
    """ATOM   1261  N   CYS A 162      -1.818  -8.192  19.088  1.00 34.25           N  
ATOM   1262  CA  CYS A 162      -2.940  -7.276  19.447  1.00 33.59           C  
ATOM   1263  C   CYS A 162      -3.006  -7.075  21.019  1.00 34.17           C  
ATOM   1264  O   CYS A 162      -4.031  -6.661  21.609  1.00 35.47           O  
ATOM   1265  CB  CYS A 162      -4.284  -7.612  18.619  1.00 33.52           C  
ATOM   1266  SG  CYS A 162      -4.609  -6.534  16.941  1.00 28.88           S  
ATOM         H   CYS A 162      -1.117  -7.774  18.510  1.00 34.25           H  
ATOM         HA  CYS A 162      -2.759  -6.249  19.096  1.00 33.59           H  
ATOM         HB2 CYS A 162      -4.261  -8.677  18.344  1.00 33.52           H  
ATOM         HB3 CYS A 162      -5.145  -7.472  19.290  1.00 33.52           H  
ATOM         HG  CYS A 162      -5.721  -6.918  16.387  1.00 28.88           H  
""", {
        Feature('HDON', [-0.9768, -7.69, 18.39]),
        Feature('HACC', [-4.685, -6.397, 21.99]),
        Feature('LIPO', [-3.282, -6.307, 15.3]),
        Feature('LIPO', [-4.244, -9.451, 18.14]),
        Feature('LIPO', [-4.674, -4.416, 16.78]),
        Feature('LIPO', [-5.771, -7.37, 19.78]),
        Feature('LIPO', [-6.405, -7.154, 16.05]),
    },
    id='CYS'
)

CYS_TRUNCATED = pytest.param(
    # Regression for CYS220 of 2ZG0, where S has no proton attached
    """ATOM   2020  N   CYS H 220       8.230 -15.806  23.838  1.00 24.28           N  
ATOM   2021  CA  CYS H 220       7.517 -14.534  23.758  1.00 20.23           C  
ATOM   2022  C   CYS H 220       6.238 -14.685  24.579  1.00 25.92           C  
ATOM   2023  O   CYS H 220       5.504 -15.657  24.378  1.00 24.17           O  
ATOM   2024  CB  CYS H 220       7.264 -14.122  22.306  1.00 18.59           C  
ATOM   2025  SG  CYS H 220       8.688 -14.134  21.195  1.00  23.1           S  
ATOM         HB2 CYS H 220       6.498 -14.792  21.888  1.00 18.59           H  
ATOM         HB3 CYS H 220       6.839 -13.108  22.307  1.00 18.59           H  
ATOM         H   CYS H 220       7.643 -16.611  23.751  1.00 24.28           H  
ATOM         HA  CYS H 220       8.119 -13.714  24.175  1.00 20.23           H  
""", {
        Feature('HACC', [5.028, -16.29, 24.25]),
        Feature('HDON', [7.526, -16.77, 23.73]),
        Feature('LIPO', [5.941, -15.28, 21.58]),
        Feature('LIPO', [6.53, -12.37, 22.31])
    },
    id='CYS_TRUNCATED'
)

HIS = pytest.param(
    """ATOM   1149  N   HIS A 148      -9.738   5.129  28.109  1.00 26.72           N  
ATOM   1150  CA  HIS A 148      -8.894   5.168  26.917  1.00  25.7           C  
ATOM   1151  C   HIS A 148      -9.023   6.505  26.193  1.00 26.61           C  
ATOM   1152  O   HIS A 148      -8.008   7.146  25.902  1.00 25.52           O  
ATOM   1153  CB  HIS A 148      -9.195   3.992  25.976  1.00 24.92           C  
ATOM   1154  CG  HIS A 148      -8.247   3.899  24.820  1.00 24.48           C  
ATOM   1155  ND1 HIS A 148      -7.335   2.865  24.677  1.00 24.63           N  
ATOM   1156  CD2 HIS A 148      -8.055   4.719  23.759  1.00 21.42           C  
ATOM   1157  CE1 HIS A 148      -6.622   3.062  23.578  1.00 21.47           C  
ATOM   1158  NE2 HIS A 148      -7.055   4.164  22.995  1.00 23.25           N  
ATOM         H   HIS A 148     -10.510   4.494  28.066  1.00 26.72           H  
ATOM         HA  HIS A 148      -7.850   5.066  27.247  1.00  25.7           H  
ATOM         HB2 HIS A 148      -9.156   3.054  26.549  1.00 24.92           H  
ATOM         HB3 HIS A 148     -10.221   4.092  25.592  1.00 24.92           H  
ATOM         HD1 HIS A 148      -7.232   2.093  25.304  1.00 24.63           H  
ATOM         HD2 HIS A 148      -8.597   5.653  23.548  1.00 21.42           H  
ATOM         HE1 HIS A 148      -5.808   2.417  23.214  1.00 21.47           H  
ATOM         HE2 HIS A 148      -6.710   4.537  22.133  1.00 23.25           H  
""",
    {
        Feature('AROM', [-10.18, 1.609, 21.95]),
        Feature('AROM', [-4.742, 5.875, 25.98]),
        Feature('HDON', [-10.66, 4.367, 28.06]),
        Feature('HDON', [-6.641, 4.612, 21.96]),
        Feature('HDON', [-7.211, 1.939, 25.43]),
        Feature('HACC', [-7.351, 7.561, 25.71]),
        Feature('POSC', [-5.586, 2.241, 23.11]),
        Feature('POSC', [-6.607, 4.649, 21.87]),
        Feature('POSC', [-7.201, 1.861, 25.49]),
        Feature('POSC', [-8.745, 5.908, 23.49]),
    },
    id='HIS'
)

GLU = pytest.param(
    """ATOM    519  N   GLU A  71       3.162   7.357  30.413  1.00 22.84           N  
ATOM    520  CA  GLU A  71       3.135   6.262  29.428  1.00  24.0           C  
ATOM    521  C   GLU A  71       3.877   5.061  30.020  1.00 24.17           C  
ATOM    522  O   GLU A  71       3.371   3.954  29.987  1.00 25.01           O  
ATOM    523  CB  GLU A  71       3.769   6.744  28.117  1.00  25.0           C  
ATOM    524  CG  GLU A  71       4.231   5.683  27.147  1.00 25.27           C  
ATOM    525  CD  GLU A  71       3.121   5.187  26.220  1.00 26.76           C  
ATOM    526  OE1 GLU A  71       3.481   4.598  25.155  1.00 27.17           O  
ATOM    527  OE2 GLU A  71       1.913   5.358  26.552  1.00 20.67           O  
ATOM         H   GLU A  71       3.719   8.142  30.143  1.00 22.84           H  
ATOM         HA  GLU A  71       2.103   5.955  29.202  1.00  24.0           H  
ATOM         HB2 GLU A  71       3.040   7.385  27.600  1.00  25.0           H  
ATOM         HB3 GLU A  71       4.633   7.377  28.368  1.00  25.0           H  
ATOM         HG2 GLU A  71       5.054   6.086  26.539  1.00 25.27           H  
ATOM         HG3 GLU A  71       4.635   4.830  27.712  1.00 25.27           H  
""", {
        Feature('NEGC', [3.868, 4.822, 23.49]), Feature('HACC', [3.708, 4.227, 24.48]),
        Feature('HACC', [1.149, 5.466, 26.76]), Feature('HACC', [3.039, 3.227, 29.97]),
        Feature('NEGC', [0.5654, 6.421, 26.43]), Feature('NEGC', [0.5855, 4.674, 27.41]),
        Feature('HDON', [3.83, 8.299, 30.09]), Feature('NEGC', [2.105, 4.733, 25.37]),
        Feature('NEGC', [3.888, 3.075, 24.47])
    },
    id='GLU'
)

GLN = pytest.param(
    """ATOM   2472  N   GLN A 325      -2.743   6.444  40.145  1.00 36.49           N  
ATOM   2473  CA  GLN A 325      -2.507   7.591  39.294  1.00 37.39           C  
ATOM   2474  C   GLN A 325      -2.654   8.950  40.032  1.00  36.9           C  
ATOM   2475  O   GLN A 325      -3.113   9.953  39.479  1.00 36.46           O  
ATOM   2476  CB  GLN A 325      -3.306   7.473  37.989  1.00 38.22           C  
ATOM   2477  CG  GLN A 325      -4.742   7.090  38.124  1.00  41.3           C  
ATOM   2478  CD  GLN A 325      -5.311   6.550  36.810  1.00 42.92           C  
ATOM   2479  OE1 GLN A 325      -5.072   7.112  35.723  1.00  43.0           O  
ATOM   2480  NE2 GLN A 325      -6.073   5.464  36.906  1.00  41.0           N  
ATOM         H   GLN A 325      -3.659   6.366  40.539  1.00 36.49           H  
ATOM         HA  GLN A 325      -1.446   7.583  39.004  1.00 37.39           H  
ATOM         HB2 GLN A 325      -3.254   8.438  37.464  1.00 38.22           H  
ATOM         HB3 GLN A 325      -2.810   6.731  37.347  1.00 38.22           H  
ATOM         HG2 GLN A 325      -4.846   6.327  38.910  1.00  41.3           H  
ATOM         HG3 GLN A 325      -5.327   7.965  38.445  1.00  41.3           H  
ATOM        HE21 GLN A 325      -6.238   5.046  37.799  1.00  41.0           H  
ATOM        HE22 GLN A 325      -6.481   5.065  36.085  1.00  41.0           H  
""",
    {
        Feature('HACC', [-4.919, 7.473, 35.03]), Feature('HDON', [-6.563, 4.985, 35.92]),
        Feature('HDON', [-6.271, 4.962, 37.98]), Feature('HACC', [-3.411, 10.6, 39.12]),
        Feature('HDON', [-3.842, 6.35, 40.62])
    },
    id='GLN'
)

GLY = pytest.param(
    """ATOM    843  N   GLY A 110       1.267  -5.826  13.366  1.00 42.11           N  
ATOM    844  CA  GLY A 110       1.283  -5.806  11.912  1.00 40.92           C  
ATOM    845  C   GLY A 110       0.106  -4.969  11.461  1.00  40.3           C  
ATOM    846  O   GLY A 110       0.251  -3.784  11.175  1.00 41.29           O  
ATOM         H   GLY A 110       0.814  -6.618  13.775  1.00 42.11           H  
ATOM         HA2 GLY A 110       2.227  -5.380  11.540  1.00 40.92           H  
ATOM         HA3 GLY A 110       1.210  -6.827  11.509  1.00 40.92           H  
""",
    {
        Feature('HACC', [0.3455, -3.012, 10.99]),
        Feature('HDON', [0.7234, -6.776, 13.86])
    },
    id='GLY'
)

ILE = pytest.param(
    """ATOM    633  N   ILE A  84      -2.143  -4.094  24.952  1.00  26.1           N  
ATOM    634  CA  ILE A  84      -0.672  -3.934  24.739  1.00 25.39           C  
ATOM    635  C   ILE A  84       0.057  -4.088  26.047  1.00 25.86           C  
ATOM    636  O   ILE A  84      -0.436  -3.650  27.096  1.00 27.21           O  
ATOM    637  CB  ILE A  84      -0.249  -2.576  24.039  1.00 25.11           C  
ATOM    638  CG1 ILE A  84       1.255  -2.563  23.706  1.00 25.01           C  
ATOM    639  CG2 ILE A  84      -0.696  -1.327  24.870  1.00 25.25           C  
ATOM    640  CD1 ILE A  84       1.708  -1.502  22.736  1.00 26.18           C  
ATOM         H   ILE A  84      -2.695  -3.410  24.474  1.00  26.1           H  
ATOM         HA  ILE A  84      -0.385  -4.731  24.037  1.00 25.39           H  
ATOM         HB  ILE A  84      -0.786  -2.512  23.081  1.00 25.11           H  
ATOM        HG12 ILE A  84       1.817  -2.438  24.644  1.00 25.01           H  
ATOM        HG13 ILE A  84       1.528  -3.547  23.296  1.00 25.01           H  
ATOM        HG21 ILE A  84      -0.384  -0.409  24.350  1.00 25.25           H  
ATOM        HG22 ILE A  84      -1.791  -1.331  24.980  1.00 25.25           H  
ATOM        HG23 ILE A  84      -0.228  -1.361  25.865  1.00 25.25           H  
ATOM        HD11 ILE A  84       2.793  -1.590  22.576  1.00 26.18           H  
ATOM        HD12 ILE A  84       1.185  -1.633  21.777  1.00 26.18           H  
ATOM        HD13 ILE A  84       1.478  -0.507  23.146  1.00 26.18           H  
""",
    {
        Feature('LIPO', [-2.587, -1.334, 25.06]), Feature('LIPO', [1.311, 0.2163, 23.44]),
        Feature('LIPO', [-0.1571, 0.2585, 23.97]), Feature('HDON', [-2.805, -3.273, 24.38]),
        Feature('LIPO', [1.726, -4.262, 23.0]), Feature('HACC', [-0.7543, -3.367, 27.77]),
        Feature('LIPO', [2.226, -2.347, 25.33]), Feature('LIPO', [-1.177, -2.465, 22.38]),
        Feature('LIPO', [3.582, -1.654, 22.46]), Feature('LIPO', [0.1123, -1.386, 26.59]),
        Feature('LIPO', [-0.9244, -0.6887, 25.29]), Feature('LIPO', [0.8047, -1.728, 21.08]),
        Feature('LIPO', [1.948, -0.9389, 22.22])
    },
    id='ILE'
)

LEU = pytest.param(
    """ATOM    547  N   LEU A  74       1.830   3.390  32.767  1.00 25.69           N  
ATOM    548  CA  LEU A  74       0.650   2.685  32.235  1.00 25.37           C  
ATOM    549  C   LEU A  74       1.068   1.368  31.564  1.00 24.76           C  
ATOM    550  O   LEU A  74       0.497   0.308  31.835  1.00 25.58           O  
ATOM    551  CB  LEU A  74      -0.170   3.588  31.272  1.00 25.28           C  
ATOM    552  CG  LEU A  74      -1.460   2.936  30.772  1.00 27.48           C  
ATOM    553  CD1 LEU A  74      -2.593   3.911  30.797  1.00 29.49           C  
ATOM    554  CD2 LEU A  74      -1.299   2.243  29.396  1.00 30.32           C  
ATOM         H   LEU A  74       2.092   4.206  32.252  1.00 25.69           H  
ATOM         HA  LEU A  74      -0.011   2.440  33.079  1.00 25.37           H  
ATOM         HB2 LEU A  74      -0.420   4.527  31.787  1.00 25.28           H  
ATOM         HB3 LEU A  74       0.457   3.852  30.407  1.00 25.28           H  
ATOM         HG  LEU A  74      -1.704   2.124  31.473  1.00 27.48           H  
ATOM        HD11 LEU A  74      -3.508   3.420  30.434  1.00 29.49           H  
ATOM        HD12 LEU A  74      -2.751   4.264  31.827  1.00 29.49           H  
ATOM        HD13 LEU A  74      -2.355   4.767  30.149  1.00 29.49           H  
ATOM        HD21 LEU A  74      -2.258   1.796  29.095  1.00 30.32           H  
ATOM        HD22 LEU A  74      -0.989   2.985  28.645  1.00 30.32           H  
ATOM        HD23 LEU A  74      -0.535   1.455  29.469  1.00 30.32           H  
""",
    {
        Feature('LIPO', [-3.199, 4.433, 30.81]), Feature('LIPO', [-2.955, 1.471, 28.88]),
        Feature('LIPO', [-0.7636, 3.524, 28.1]), Feature('LIPO', [-0.6019, 5.21, 32.16]),
        Feature('LIPO', [-2.182, 5.39, 29.68]), Feature('LIPO', [0.02064, 0.8819, 29.52]),
        Feature('LIPO', [-2.866, 4.521, 32.58]), Feature('HACC', [0.1269, -0.3791, 32.01]),
        Feature('HDON', [2.144, 4.369, 32.15]), Feature('LIPO', [-4.173, 3.063, 30.17]),
        Feature('LIPO', [-1.216, 1.885, 28.69]), Feature('LIPO', [0.9128, 4.044, 29.78]),
        Feature('LIPO', [-1.881, 1.534, 31.98])
    },
    id='LEU'
)

LYS = pytest.param(
    """ATOM    362  N   LYS A  53      10.305   1.187  20.473  1.00 33.62           N  
ATOM    363  CA  LYS A  53       9.944   2.582  20.279  1.00 33.93           C  
ATOM    364  C   LYS A  53      10.939   3.457  21.043  1.00 33.03           C  
ATOM    365  O   LYS A  53      11.151   3.278  22.245  1.00  32.6           O  
ATOM    366  CB  LYS A  53       8.487   2.843  20.699  1.00 33.99           C  
ATOM    367  CG  LYS A  53       8.015   4.329  20.683  1.00 36.22           C  
ATOM    368  CD  LYS A  53       6.534   4.444  21.124  1.00 35.43           C  
ATOM    369  CE  LYS A  53       6.114   5.839  21.704  1.00 38.16           C  
ATOM    370  NZ  LYS A  53       5.339   5.781  23.107  1.00 32.75           N  
ATOM         H   LYS A  53      10.159   0.841  21.400  1.00 33.62           H  
ATOM         HA  LYS A  53      10.001   2.839  19.211  1.00 33.93           H  
ATOM         HB2 LYS A  53       7.828   2.265  20.034  1.00 33.99           H  
ATOM         HB3 LYS A  53       8.345   2.447  21.716  1.00 33.99           H  
ATOM         HG2 LYS A  53       8.649   4.927  21.354  1.00 36.22           H  
ATOM         HG3 LYS A  53       8.134   4.745  19.672  1.00 36.22           H  
ATOM         HD2 LYS A  53       5.893   4.217  20.260  1.00 35.43           H  
ATOM         HD3 LYS A  53       6.334   3.674  21.884  1.00 35.43           H  
ATOM         HE2 LYS A  53       7.017   6.456  21.826  1.00 38.16           H  
ATOM         HE3 LYS A  53       5.473   6.348  20.969  1.00 38.16           H  
ATOM         HZ1 LYS A  53       4.820   6.626  23.237  1.00 32.75           H  
ATOM         HZ2 LYS A  53       4.712   5.002  23.110  1.00 32.75           H  
ATOM         HZ3 LYS A  53       6.001   5.679  23.850  1.00 32.75           H  
""",
    {
        Feature('POSC', [4.524, 4.768, 23.11]), Feature('LIPO', [5.427, 4.052, 19.63]),
        Feature('LIPO', [6.189, 3.114, 22.44]), Feature('LIPO', [8.221, 5.048, 18.94]),
        Feature('POSC', [4.664, 6.879, 23.28]), Feature('LIPO', [7.349, 1.845, 19.55]),
        Feature('POSC', [5.194, 5.77, 23.37]), Feature('LIPO', [8.242, 2.159, 22.46]),
        Feature('HDON', [10.13, 0.7718, 21.59]), Feature('HDON', [4.716, 6.795, 23.26]),
        Feature('LIPO', [9.11, 5.362, 21.84]), Feature('HDON', [4.587, 4.846, 23.11]),
        Feature('POSC', [6.2, 5.648, 24.07]), Feature('HDON', [6.133, 5.659, 24.0]),
        Feature('HACC', [11.29, 3.162, 23.02])
    },
    id='LYS'
)

MET = pytest.param(
    """ATOM    835  N   MET A 109       2.893  -4.854  16.219  1.00 43.58           N  
ATOM    836  CA  MET A 109       1.578  -4.786  15.550  1.00  43.5           C  
ATOM    837  C   MET A 109       1.828  -4.839  14.051  1.00 42.49           C  
ATOM    838  O   MET A 109       2.518  -3.998  13.515  1.00  42.6           O  
ATOM    839  CB  MET A 109       0.770  -3.487  15.882  1.00 43.44           C  
ATOM    840  CG  MET A 109      -0.429  -3.188  14.862  1.00 43.83           C  
ATOM    841  SD  MET A 109      -1.534  -1.676  14.815  1.00 45.32           S  
ATOM    842  CE  MET A 109      -0.542  -0.537  13.837  1.00 45.08           C  
ATOM         H   MET A 109       3.247  -3.958  16.488  1.00 43.58           H  
ATOM         HA  MET A 109       0.971  -5.630  15.911  1.00  43.5           H  
ATOM         HB2 MET A 109       0.360  -3.572  16.899  1.00 43.44           H  
ATOM         HB3 MET A 109       1.459  -2.629  15.884  1.00 43.44           H  
ATOM         HG2 MET A 109       0.021  -3.249  13.860  1.00 43.83           H  
ATOM         HG3 MET A 109      -1.113  -4.044  14.960  1.00 43.83           H  
ATOM         HE1 MET A 109      -1.001   0.462  13.860  1.00 45.08           H  
ATOM         HE2 MET A 109       0.474  -0.483  14.255  1.00 45.08           H  
ATOM         HE3 MET A 109      -0.492  -0.893  12.797  1.00 45.08           H  
""",
    {
        Feature('LIPO', [1.213, -0.4437, 14.56]), Feature('HACC', [2.973, -3.443, 13.16]),
        Feature('LIPO', [1.96, -2.005, 15.89]), Feature('LIPO', [-0.4556, -1.152, 12.04]),
        Feature('LIPO', [-2.581, -0.9279, 16.66]), Feature('LIPO', [-3.741, -1.845, 14.42]),
        Feature('LIPO', [-1.335, 1.189, 13.88]), Feature('LIPO', [0.3482, -3.293, 13.13]),
        Feature('HDON', [3.318, -3.779, 16.54]), Feature('LIPO', [-1.61, -4.666, 15.03]),
        Feature('LIPO', [0.06177, -3.634, 17.64]), Feature('LIPO', [-0.101, -0.03061, 13.4]),
    },
    id='MET'
)

PHE = pytest.param(
    """ATOM   1317  N   PHE A 169      -1.700   4.733  20.035  1.00 34.89           N  
ATOM   1318  CA  PHE A 169      -0.694   5.577  19.401  1.00 38.89           C  
ATOM   1319  C   PHE A 169      -0.707   7.078  19.808  1.00 40.29           C  
ATOM   1320  O   PHE A 169      -0.108   7.461  20.802  1.00 42.12           O  
ATOM   1321  CB  PHE A 169      -0.866   5.456  17.886  1.00 38.54           C  
ATOM   1322  CG  PHE A 169      -0.443   4.147  17.343  1.00 39.88           C  
ATOM   1323  CD1 PHE A 169       0.863   3.713  17.493  1.00 38.52           C  
ATOM   1324  CD2 PHE A 169      -1.336   3.357  16.642  1.00 40.94           C  
ATOM   1325  CE1 PHE A 169       1.267   2.518  16.973  1.00 39.97           C  
ATOM   1326  CE2 PHE A 169      -0.937   2.158  16.110  1.00 40.63           C  
ATOM   1327  CZ  PHE A 169       0.364   1.728  16.280  1.00 41.45           C  
ATOM         H   PHE A 169      -2.390   4.382  19.402  1.00 34.89           H  
ATOM         HA  PHE A 169       0.280   5.208  19.754  1.00 38.89           H  
ATOM         HB2 PHE A 169      -1.923   5.625  17.631  1.00 38.54           H  
ATOM         HB3 PHE A 169      -0.286   6.251  17.395  1.00 38.54           H  
ATOM         HD1 PHE A 169       1.584   4.340  18.038  1.00 38.52           H  
ATOM         HD2 PHE A 169      -2.375   3.693  16.510  1.00 40.94           H  
ATOM         HE1 PHE A 169       2.307   2.183  17.104  1.00 39.97           H  
ATOM         HE2 PHE A 169      -1.652   1.539  15.548  1.00 40.63           H  
ATOM         HZ  PHE A 169       0.683   0.760  15.866  1.00 41.45           H  
""",
    {
        Feature('AROM', [-1.038, 1.104, 20.22]),
        Feature('AROM', [0.964, 4.77, 13.4]),
        Feature('HDON', [-2.528, 4.312, 19.28]),
        Feature('HACC', [0.2841, 7.712, 21.45]),
        Feature('LIPO', [-0.2372, 2.57, 17.49]),
        Feature('LIPO', [-2.172, 1.089, 15.14]),
        Feature('LIPO', [-3.131, 3.937, 16.41]),
        Feature('LIPO', [0.1632, 3.303, 16.12]),
        Feature('LIPO', [0.915, 0.05605, 15.56]),
        Feature('LIPO', [2.108, 4.796, 18.43]),
        Feature('LIPO', [3.063, 1.939, 17.2]),
    },
    id='PHE'
)

PRO = pytest.param(
    """ATOM    218  N   PRO A  29      13.094   0.610  10.492  1.00 51.86           N  
ATOM    219  CA  PRO A  29      11.999   1.586  10.589  1.00 52.19           C  
ATOM    220  C   PRO A  29      10.605   1.049  10.226  1.00 52.44           C  
ATOM    221  O   PRO A  29      10.462   0.326   9.247  1.00 52.97           O  
ATOM    222  CB  PRO A  29      12.448   2.681   9.626  1.00 52.37           C  
ATOM    223  CG  PRO A  29      13.933   2.686   9.773  1.00 51.81           C  
ATOM    224  CD  PRO A  29      14.347   1.257  10.044  1.00 51.47           C  
ATOM         HA  PRO A  29      11.852   1.915  11.628  1.00 52.19           H  
ATOM         HB2 PRO A  29      12.144   2.461   8.592  1.00 52.37           H  
ATOM         HB3 PRO A  29      12.012   3.656   9.889  1.00 52.37           H  
ATOM         HG2 PRO A  29      14.415   3.066   8.860  1.00 51.81           H  
ATOM         HG3 PRO A  29      14.242   3.345  10.598  1.00 51.81           H  
ATOM         HD2 PRO A  29      14.751   0.775   9.142  1.00 51.47           H  
ATOM         HD3 PRO A  29      15.129   1.202  10.816  1.00 51.47           H  
""",
    {
        Feature('LIPO', [14.47, 3.824, 11.2]), Feature('HACC', [10.37, -0.146, 8.608]),
        Feature('LIPO', [14.77, 3.342, 8.196]), Feature('LIPO', [11.69, 4.365, 10.08]),
        Feature('LIPO', [15.7, 1.162, 11.38]), Feature('LIPO', [11.92, 2.301, 7.84]),
        Feature('LIPO', [15.04, 0.4243, 8.486])
    },
    id='PRO'
)

SER = pytest.param(
    """ATOM    212  N   SER A  28      14.655  -2.349  11.767  1.00  50.2           N  
ATOM    213  CA  SER A  28      14.214  -1.583  10.590  1.00 50.93           C  
ATOM    214  C   SER A  28      12.992  -0.688  10.849  1.00 51.08           C  
ATOM    215  O   SER A  28      11.997  -1.142  11.398  1.00 51.17           O  
ATOM    216  CB  SER A  28      14.012  -2.480   9.351  1.00 51.07           C  
ATOM    217  OG  SER A  28      12.838  -3.256   9.437  1.00 52.66           O  
ATOM         H   SER A  28      13.913  -2.730  12.319  1.00  50.2           H  
ATOM         HA  SER A  28      15.044  -0.895  10.371  1.00 50.93           H  
ATOM         HB2 SER A  28      13.968  -1.851   8.449  1.00 51.07           H  
ATOM         HB3 SER A  28      14.881  -3.145   9.237  1.00 51.07           H  
ATOM         HG  SER A  28      12.748  -3.818   8.615  1.00 52.66           H  
""",
    {
        Feature('HDON', [13.76, -2.806, 12.43]), Feature('HACC', [11.35, -1.439, 11.76]),
        Feature('HACC', [12.74, -3.85, 9.963]), Feature('HDON', [12.73, -3.93, 8.451]),
        Feature('HACC', [12.1, -2.948, 9.416])
    },
    id='SER'
)

THR = pytest.param(
    """ATOM    810  N   THR A 106       7.772  -3.530  22.865  1.00 34.38           N  
ATOM    811  CA  THR A 106       6.393  -3.929  22.772  1.00 38.11           C  
ATOM    812  C   THR A 106       6.361  -5.072  21.725  1.00  39.7           C  
ATOM    813  O   THR A 106       7.324  -5.239  20.966  1.00 40.19           O  
ATOM    814  CB  THR A 106       5.600  -2.683  22.302  1.00 38.17           C  
ATOM    815  OG1 THR A 106       4.969  -2.034  23.424  1.00 39.56           O  
ATOM    816  CG2 THR A 106       4.611  -3.011  21.253  1.00 37.97           C  
ATOM         H   THR A 106       8.167  -3.291  21.978  1.00 34.38           H  
ATOM         HA  THR A 106       5.953  -4.284  23.715  1.00 38.11           H  
ATOM         HB  THR A 106       6.321  -1.985  21.851  1.00 38.17           H  
ATOM         HG1 THR A 106       4.896  -1.053  23.244  1.00 39.56           H  
ATOM        HG21 THR A 106       4.076  -2.097  20.954  1.00 37.97           H  
ATOM        HG22 THR A 106       5.129  -3.435  20.380  1.00 37.97           H  
ATOM        HG23 THR A 106       3.891  -3.745  21.644  1.00 37.97           H  
""",
    {
        Feature('HDON', [4.881, -0.8568, 23.21]),
        Feature('HDON', [8.246, -3.243, 21.8]),
        Feature('HACC', [4.204, -2.134, 23.63]),
        Feature('HACC', [5.299, -1.958, 24.15]),
        Feature('HACC', [7.947, -5.347, 20.48]),
        Feature('LIPO', [3.367, -4.279, 21.93]),
        Feature('LIPO', [3.687, -1.433, 20.74]),
        Feature('LIPO', [4.076, -3.188, 20.69]),
        Feature('LIPO', [5.506, -3.743, 19.75]),
    },
    id='THR'
)

TRP = pytest.param(
    """ATOM    118  N   TRP A  18      20.313   5.156  17.273  1.00 52.47           N  
ATOM    119  CA  TRP A  18      19.704   4.118  18.098  1.00  50.2           C  
ATOM    120  C   TRP A  18      20.747   3.390  18.911  1.00 49.86           C  
ATOM    121  O   TRP A  18      21.686   2.834  18.348  1.00 50.21           O  
ATOM    122  CB  TRP A  18      18.936   3.119  17.237  1.00 48.82           C  
ATOM    123  CG  TRP A  18      17.615   3.658  16.755  1.00 47.99           C  
ATOM    124  CD1 TRP A  18      17.429   4.644  15.828  1.00 46.59           C  
ATOM    125  CD2 TRP A  18      16.297   3.244  17.168  1.00 45.48           C  
ATOM    126  NE1 TRP A  18      16.089   4.883  15.651  1.00 46.99           N  
ATOM    127  CE2 TRP A  18      15.368   4.035  16.449  1.00 46.37           C  
ATOM    128  CE3 TRP A  18      15.815   2.298  18.075  1.00 44.64           C  
ATOM    129  CZ2 TRP A  18      13.978   3.903  16.607  1.00 46.37           C  
ATOM    130  CZ3 TRP A  18      14.428   2.161  18.226  1.00 46.66           C  
ATOM    131  CH2 TRP A  18      13.533   2.959  17.495  1.00 46.51           C  
ATOM         H   TRP A  18      20.777   4.809  16.458  1.00 52.47           H  
ATOM         HA  TRP A  18      19.001   4.612  18.785  1.00  50.2           H  
ATOM         HB2 TRP A  18      19.552   2.841  16.369  1.00 48.82           H  
ATOM         HB3 TRP A  18      18.762   2.200  17.816  1.00 48.82           H  
ATOM         HD1 TRP A  18      18.237   5.171  15.299  1.00 46.59           H  
ATOM         HE1 TRP A  18      15.700   5.569  15.036  1.00 46.99           H  
ATOM         HE3 TRP A  18      16.508   1.674  18.658  1.00 44.64           H  
ATOM         HZ2 TRP A  18      13.275   4.532  16.041  1.00 46.37           H  
ATOM         HZ3 TRP A  18      14.033   1.414  18.931  1.00 46.66           H  
ATOM         HH2 TRP A  18      12.451   2.823  17.638  1.00 46.51           H  
""",
    {
        Feature('AROM', [14.38, 0.1145, 14.66]),
        Feature('AROM', [14.49, 5.546, 20.54]),
        Feature('AROM', [16.5, 1.377, 13.43]),
        Feature('AROM', [16.62, 6.809, 19.31]),
        Feature('HDON', [15.62, 5.706, 14.91]),
        Feature('HDON', [20.87, 4.74, 16.29]),
        Feature('HACC', [22.3, 2.472, 17.98]),
        Feature('LIPO', [11.66, 2.724, 17.74]),
        Feature('LIPO', [12.76, 4.989, 15.63]),
        Feature('LIPO', [13.75, 0.871, 19.44]),
        Feature('LIPO', [14.43, 2.287, 17.01]),
        Feature('LIPO', [14.45, 3.373, 18.19]),
        Feature('LIPO', [17.01, 1.22, 19.08]),
    },
    id='TRP'
)

TYR = pytest.param(
    """ATOM   1084  N   TYR A 140     -10.867  -4.696  31.686  1.00  28.5           N  
ATOM   1085  CA  TYR A 140      -9.622  -4.665  32.449  1.00 26.65           C  
ATOM   1086  C   TYR A 140      -9.137  -3.234  32.681  1.00 27.02           C  
ATOM   1087  O   TYR A 140      -8.759  -2.881  33.800  1.00 28.11           O  
ATOM   1088  CB  TYR A 140      -8.573  -5.543  31.753  1.00 25.69           C  
ATOM   1089  CG  TYR A 140      -7.166  -5.472  32.334  1.00 24.08           C  
ATOM   1090  CD1 TYR A 140      -6.793  -6.269  33.418  1.00  24.5           C  
ATOM   1091  CD2 TYR A 140      -6.205  -4.591  31.783  1.00 22.23           C  
ATOM   1092  CE1 TYR A 140      -5.490  -6.224  33.938  1.00 24.92           C  
ATOM   1093  CE2 TYR A 140      -4.919  -4.516  32.293  1.00 21.42           C  
ATOM   1094  CZ  TYR A 140      -4.572  -5.318  33.386  1.00 24.72           C  
ATOM   1095  OH  TYR A 140      -3.303  -5.232  33.917  1.00 27.93           O  
ATOM         H   TYR A 140     -10.794  -5.171  30.809  1.00  28.5           H  
ATOM         HA  TYR A 140      -9.802  -5.082  33.451  1.00 26.65           H  
ATOM         HB2 TYR A 140      -8.913  -6.588  31.790  1.00 25.69           H  
ATOM         HB3 TYR A 140      -8.527  -5.257  30.692  1.00 25.69           H  
ATOM         HD1 TYR A 140      -7.533  -6.944  33.872  1.00  24.5           H  
ATOM         HD2 TYR A 140      -6.485  -3.953  30.932  1.00 22.23           H  
ATOM         HE1 TYR A 140      -5.194  -6.887  34.764  1.00 24.92           H  
ATOM         HE2 TYR A 140      -4.181  -3.835  31.845  1.00 21.42           H  
ATOM         HH  TYR A 140      -3.345  -5.374  34.906  1.00 27.93           H  
""",
    {
        Feature('AROM', [-4.644, -8.202, 30.28]),
        Feature('AROM', [-7.071, -2.595, 35.44]),
        Feature('HDON', [-10.78, -5.266, 30.63]),
        Feature('HDON', [-3.353, -5.402, 35.1]),
        Feature('HACC', [-2.732, -5.77, 33.76]),
        Feature('HACC', [-2.885, -4.55, 33.93]),
        Feature('HACC', [-8.514, -2.652, 34.53]),
    },
    id='TYR'
)

# 3HEG has a TYR on A 35 which is truncated
TYR_TRUNCATED = pytest.param(
    """ATOM    232  N   TYR A  35       3.890   5.530  15.364  1.00 56.36           N  
ATOM    233  CA  TYR A  35       4.070   6.446  16.490  1.00 55.91           C  
ATOM    234  C   TYR A  35       5.048   7.601  16.179  1.00  55.3           C  
ATOM    235  O   TYR A  35       4.746   8.521  15.397  1.00 55.26           O  
ATOM    236  CB  TYR A  35       4.508   5.658  17.762  1.00 55.74           C  
ATOM         HA  TYR A  35       3.096   6.919  16.681  1.00 55.91           H  
""",
    {
        Feature('HACC', [4.552, 9.112, 14.89])
    },
    id='TYR_TRUNCATED'
)

VAL = pytest.param(
    """ATOM    247  N   VAL A  38       9.777   4.049  15.622  1.00 43.74           N  
ATOM    248  CA  VAL A  38       9.605   2.654  15.941  1.00 42.34           C  
ATOM    249  C   VAL A  38      10.431   1.869  14.930  1.00 42.26           C  
ATOM    250  O   VAL A  38      10.499   2.213  13.730  1.00 40.89           O  
ATOM    251  CB  VAL A  38       8.090   2.222  15.901  1.00 42.72           C  
ATOM    252  CG1 VAL A  38       7.916   0.747  16.163  1.00 40.63           C  
ATOM    253  CG2 VAL A  38       7.263   3.017  16.922  1.00 41.47           C  
ATOM         H   VAL A  38       9.553   4.267  14.672  1.00 43.74           H  
ATOM         HA  VAL A  38       9.943   2.453  16.968  1.00 42.34           H  
ATOM         HB  VAL A  38       7.729   2.441  14.885  1.00 42.72           H  
ATOM        HG11 VAL A  38       6.847   0.491  16.126  1.00 40.63           H  
ATOM        HG12 VAL A  38       8.457   0.172  15.397  1.00 40.63           H  
ATOM        HG13 VAL A  38       8.317   0.501  17.157  1.00 40.63           H  
ATOM        HG21 VAL A  38       6.212   2.696  16.873  1.00 41.47           H  
ATOM        HG22 VAL A  38       7.655   2.835  17.934  1.00 41.47           H  
ATOM        HG23 VAL A  38       7.329   4.091  16.692  1.00 41.47           H  
    """,
    {
        Feature('LIPO', [6.832, 3.431, 17.45]), Feature('LIPO', [6.069, 0.3048, 16.1]),
        Feature('HACC', [10.54, 2.433, 12.96]), Feature('LIPO', [7.824, -0.03541, 16.3]),
        Feature('LIPO', [7.467, 2.6, 14.15]), Feature('LIPO', [7.377, 4.872, 16.52]),
        Feature('LIPO', [8.85, -0.2462, 14.84]), Feature('LIPO', [8.609, 0.322, 17.88]),
        Feature('LIPO', [7.94, 2.703, 18.67]), Feature('LIPO', [5.448, 2.463, 16.84]),
        Feature('HDON', [9.508, 4.311, 14.48])
    },
    id='VAL'
)


def feat2point(f: Feature):
    return f.kind, f.position[0], f.position[1], f.position[2]


def feat2points(features):
    return [feat2point(f) for f in features]


residues = [
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    CYS_TRUNCATED,
    HIS,
    GLU,
    GLN,
    GLY,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    TYR_TRUNCATED,
    VAL,
]


@pytest.mark.parametrize("block,expected", residues)
def test_from_site(block, expected, ligand_3heg_bax: Ligand):
    site = prep_site(block, ligand_3heg_bax.molecule)

    features = from_site(site)

    assert_features(expected, features)


def dump_artificial_pdb_with_all_residues(all_res_site_pdb='all_res_site.pdb'):
    """

    ```python
    from tests.pharmacohore.test_from_site import dump_artificial_pdb_with_all_residues
    dump_artificial_pdb_with_all_residues('all_res_site_pdb')
    ```

    """
    ligand = ligands(pdb_from_file('tests/fixtures/3HEG.prepped.pdb'))[0]
    block = ""
    for residue in residues:
        block += residue.values[0]
    site = prep_site(block, ligand.molecule)
    block_with_ligand = site.ligand().model().to_file_string('pdb')
    with open(all_res_site_pdb, 'w') as f:
        f.write(block_with_ligand)


def dump_pharmacophore_with_all_residues(all_res_site_phar='all_res_site.phar'):
    """"

    ```python
    from tests.pharmacohore.test_from_site import dump_pharmacophore_with_all_residues
    dump_pharmacophore_with_all_residues('all_res_site.phar')
    ```

    """
    features = set()
    for residue in residues:
        features |= residue.values[1]
    phar = as_phar(all_res_site_phar, feat2points(features))
    with open(all_res_site_phar, 'w') as f:
        f.write(phar)


def dump4molviewer(filename='kripo-phar-molviewer.json'):
    """Dumps each pharmacophore/residue pairs to a json file which can be viewed with the 3d-e-chem molviewer

    ```python
    from tests.pharmacohore.test_from_site import dump4molviewer
    dump4molviewer()
    ```

    """
    rows = []
    ligand = ligands(pdb_from_file('tests/fixtures/3HEG.prepped.pdb'))[0]
    for residue in residues:
        block = residue.values[0]
        features = residue.values[1]
        site = prep_site(block, ligand.molecule)
        label = site.residue().name() + str(len(list(site.residue().atoms())))
        row = {
            'id': label,
            'label': label,
            'protein': {
                'data': site.ligand().model().to_file_string('pdb'),
                'format': 'pdb',
            },
            'pharmacophore': {
                'data': as_phar(label, feat2points(features)),
                'format': 'phar',
            },
        }
        rows.append(row)
    with open(filename, 'w') as f:
        json.dump(rows, f)
