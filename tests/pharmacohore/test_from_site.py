from typing import Set
import json

from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.structures.chains import Site
import pytest

from kripo.ligand import Ligand
from kripo.pharmacophore import Feature, from_site
from kripodb.pharmacophores import as_phar


def prep_site(block, ligand):
    model = pdb_dict_to_pdb(pdb_string_to_pdb_dict(block)).model()
    model.add_molecule(ligand)
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
        Feature('HACC', [10.82, -1.872, 12.71]),
        Feature('HDON', [13.77, -3.462, 13.6]),
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
        Feature('HACC', [-0.1856, 14.66, 32.12]), Feature('POSC', [0.265, 12.97, 29.97]),
        Feature('POSC', [-0.921, 11.09, 31.74]), Feature('LIPO', [1.356, 10.93, 33.31]),
        Feature('HACC', [4.044, 10.14, 32.08]), Feature('HDON', [1.42, 5.654, 32.1]),
        Feature('LIPO', [-0.7883, 8.522, 33.04]), Feature('HACC', [0.2534, 12.05, 29.39]),
        Feature('LIPO', [2.389, 10.49, 30.47]), Feature('LIPO', [0.3552, 7.631, 30.35]),
        Feature('HACC', [-1.369, 11.18, 32.73]), Feature('HACC', [-1.081, 13.27, 33.38]),
        Feature('POSC', [-0.5965, 13.75, 32.54]), Feature('HACC', [0.5198, 14.02, 30.01])}
    ,
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
        Feature('HACC', [-3.649, 3.007, 18.54]), Feature('HDON', [-6.384, 5.919, 19.27]),
        Feature('HACC', [-3.511, 5.078, 18.65]), Feature('HACC', [-7.558, 2.352, 16.62]),
        Feature('HDON', [-5.815, -0.09534, 20.33])
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
        Feature('HDON', [-4.368, 7.442, 24.32]), Feature('NEGC', [-4.883, 6.607, 20.21]),
        Feature('HDON', [-0.03166, 4.629, 22.5]), Feature('HACC', [-1.68, 2.363, 22.7]),
        Feature('HDON', [-3.917, 6.769, 20.85])
    },
    id='ASP'
)

CYS = pytest.param(
    """ATOM   1261  N   CYS A 162      -6.818 -13.192  14.088  1.00 34.25           N  
ATOM   1262  CA  CYS A 162      -7.940 -12.276  14.447  1.00 33.59           C  
ATOM   1263  C   CYS A 162      -8.006 -12.075  16.019  1.00 34.17           C  
ATOM   1264  O   CYS A 162      -9.031 -11.661  16.609  1.00 35.47           O  
ATOM   1265  CB  CYS A 162      -9.284 -12.612  13.619  1.00 33.52           C  
ATOM   1266  SG  CYS A 162      -9.609 -11.534  11.941  1.00 28.88           S  
ATOM         H   CYS A 162      -6.117 -12.774  13.510  1.00 34.25           H  
ATOM         HA  CYS A 162      -7.759 -11.249  14.096  1.00 33.59           H  
ATOM         HB2 CYS A 162      -9.261 -13.677  13.344  1.00 33.52           H  
ATOM         HB3 CYS A 162     -10.145 -12.472  14.290  1.00 33.52           H  
ATOM         HG  CYS A 162     -10.721 -11.918  11.387  1.00 28.88           H  
""",
    {
        Feature('HACC', [-5.977, -12.69, 13.39]), Feature('HDON', [-11.41, -12.15, 11.05]),
        Feature('HDON', [-9.685, -11.4, 16.99]), Feature('HDON', [-9.738, -11.11, 11.28]),
        Feature('HDON', [-9.244, -14.45, 13.14]), Feature('HDON', [-10.77, -12.37, 14.78])
    },
    id='CYS'
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
        Feature('HDON', [-7.351, 7.561, 25.71]), Feature('AROM', [-9.775, 1.929, 22.26]),
        Feature('HACC', [-10.66, 4.367, 28.06]), Feature('POSC', [-6.607, 4.649, 21.87]),
        Feature('HACC', [-6.641, 4.612, 21.96]), Feature('POSC', [-7.201, 1.861, 25.49]),
        Feature('POSC', [-8.745, 5.908, 23.49]), Feature('POSC', [-5.586, 2.241, 23.11]),
        Feature('HACC', [-7.211, 1.939, 25.43]), Feature('AROM', [-5.15, 5.555, 25.68])
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
""",{
        Feature('NEGC', [3.868, 4.822, 23.49]), Feature('HDON', [3.708, 4.227, 24.48]),
        Feature('HDON', [1.149, 5.466, 26.76]), Feature('HDON', [3.039, 3.227, 29.97]),
        Feature('NEGC', [0.5654, 6.421, 26.43]), Feature('NEGC', [0.5855, 4.674, 27.41]),
        Feature('HACC', [3.83, 8.299, 30.09]), Feature('NEGC', [2.105, 4.733, 25.37]),
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
        Feature('HDON', [-4.919, 7.473, 35.03]), Feature('HACC', [-6.563, 4.985, 35.92]),
        Feature('HACC', [-6.271, 4.962, 37.98]), Feature('HDON', [-3.411, 10.6, 39.12]),
        Feature('HACC', [-3.842, 6.35, 40.62])
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
        Feature('HDON', [0.3455, -3.012, 10.99]),
        Feature('HACC', [0.7234, -6.776, 13.86])
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
        Feature('LIPO', [-0.1571, 0.2585, 23.97]), Feature('HACC', [-2.805, -3.273, 24.38]),
        Feature('LIPO', [1.726, -4.262, 23.0]), Feature('HDON', [-0.7543, -3.367, 27.77]),
        Feature('LIPO', [2.226, -2.347, 25.33]), Feature('LIPO', [-1.177, -2.465, 22.38]),
        Feature('LIPO', [3.582, -1.654, 22.46]), Feature('LIPO', [0.1123, -1.386, 26.59]),
        Feature('LIPO', [-0.9244, -0.6887, 25.29]), Feature('LIPO', [0.8047, -1.728, 21.08]),
        Feature('LIPO', [1.948, -0.9389, 22.22])
    },
    id='ILE'
)

LEU = pytest.param(
    """ATOM   2602  N   LEU A 340      11.814  10.073  38.518  1.00 30.40           N  
ATOM   2603  CA  LEU A 340      10.725   9.132  38.690  1.00 30.79           C  
ATOM   2604  C   LEU A 340      10.583   8.167  37.520  1.00 30.62           C  
ATOM   2605  O   LEU A 340      10.274   6.992  37.721  1.00 32.64           O  
ATOM   2606  CB  LEU A 340       9.435   9.910  38.866  1.00 30.61           C  
ATOM   2607  CG  LEU A 340       8.637  10.105  40.169  1.00 33.91           C  
ATOM   2608  CD1 LEU A 340       9.363   9.912  41.543  1.00 33.34           C  
ATOM   2609  CD2 LEU A 340       7.854  11.454  40.103  1.00 30.62           C  
""",
    {
        Feature('HACC', [11.94, 10.45, 37.36]),
        Feature('HDON', [10.06, 6.232, 37.87]),
        Feature('LIPO', [10.29, 11.92, 40.32]),
        Feature('LIPO', [10.92, 10.09, 42.31]),
        Feature('LIPO', [7.017, 11.72, 40.01]),
        Feature('LIPO', [7.747, 10.22, 42.63]),
        Feature('LIPO', [8.02, 8.612, 38.47]),
        Feature('LIPO', [8.867, 8.505, 38.22]),
        Feature('LIPO', [9.045, 7.902, 41.15]),
        Feature('LIPO', [9.246, 11.29, 37.65]),
        Feature('LIPO', [9.366, 9.399, 42.23]),
        Feature('LIPO', [9.826, 10.5, 37.58]),
    },
    id='LEU'
)

LYS = pytest.param(
    """ATOM   2587  N   LYS A 338      11.629  13.271  35.523  1.00 28.97           N  
ATOM   2588  CA  LYS A 338      12.685  12.487  34.904  1.00 29.96           C  
ATOM   2589  C   LYS A 338      13.301  11.445  35.866  1.00 30.29           C  
ATOM   2590  O   LYS A 338      13.566  10.323  35.437  1.00 31.18           O  
ATOM   2591  CB  LYS A 338      13.762  13.380  34.304  1.00 29.03           C  
ATOM   2592  CG  LYS A 338      14.889  12.611  33.642  1.00 31.14           C  
ATOM   2593  CD  LYS A 338      16.223  13.370  33.666  1.00 37.81           C  
ATOM   2594  CE  LYS A 338      16.521  13.999  32.320  1.00 39.10           C  
ATOM   2595  NZ  LYS A 338      17.999  14.105  32.112  1.00 38.58           N  
""",
    {
        Feature('HACC', [10.88, 13.87, 34.97]),
        Feature('HACC', [17.43, 13.41, 30.93]),
        Feature('HACC', [18.41, 15.09, 31.66]),
        Feature('HDON', [14.31, 10.3, 36.48]),
        Feature('LIPO', [13.42, 14.62, 33.04]),
        Feature('LIPO', [14.02, 14.57, 36.0]),
        Feature('LIPO', [14.36, 13.41, 31.63]),
        Feature('LIPO', [15.11, 10.9, 34.3]),
        Feature('LIPO', [16.46, 13.3, 35.53]),
        Feature('LIPO', [16.73, 11.58, 32.14]),
        Feature('POSC', [17.39, 13.37, 30.85]),
        Feature('POSC', [18.13, 13.94, 31.7]),
        Feature('POSC', [18.45, 15.18, 31.64]),
    },
    id='LYS'
)

MET = pytest.param(
    """ATOM   2027  N   MET A 268     -34.186   6.974  15.854  1.00 67.59           N  
ATOM   2028  CA  MET A 268     -34.482   6.029  14.776  1.00 66.64           C  
ATOM   2029  C   MET A 268     -34.605   4.587  15.269  1.00 65.60           C  
ATOM   2030  O   MET A 268     -34.162   4.249  16.370  1.00 65.52           O  
ATOM   2031  CB  MET A 268     -33.450   6.144  13.653  1.00 66.75           C  
ATOM   2032  CG  MET A 268     -33.572   7.441  12.860  1.00 67.06           C  
ATOM   2033  SD  MET A 268     -32.305   7.606  11.590  1.00 68.00           S  
ATOM   2034  CE  MET A 268     -32.644   9.249  10.944  1.00 68.60           C  
""",
    {
        Feature('HACC', [-34.07, 8.058, 15.55]),
        Feature('HDON', [-35.53, 3.948, 16.54]),
        Feature('LIPO', [-31.28, 10.31, 10.83]),
        Feature('LIPO', [-32.01, 8.184, 14.36]),
        Feature('LIPO', [-33.1, 10.03, 9.343]),
        Feature('LIPO', [-33.59, 10.67, 11.56]),
        Feature('LIPO', [-34.35, 8.448, 13.95]),
    },
    id='MET'
)

PHE = pytest.param(
    """ATOM   2667  N   PHE A 348      11.898  -2.176  36.616  1.00 42.77           N  
ATOM   2668  CA  PHE A 348      11.170  -3.005  35.654  1.00 43.16           C  
ATOM   2669  C   PHE A 348      11.523  -4.480  35.867  1.00 44.16           C  
ATOM   2670  O   PHE A 348      12.692  -4.843  35.976  1.00 44.07           O  
ATOM   2671  CB  PHE A 348      11.426  -2.573  34.196  1.00 41.65           C  
ATOM   2672  CG  PHE A 348      10.790  -3.486  33.162  1.00 41.53           C  
ATOM   2673  CD1 PHE A 348       9.414  -3.518  32.989  1.00 39.20           C  
ATOM   2674  CD2 PHE A 348      11.578  -4.313  32.363  1.00 40.12           C  
ATOM   2675  CE1 PHE A 348       8.838  -4.364  32.061  1.00 38.12           C  
ATOM   2676  CE2 PHE A 348      11.006  -5.163  31.427  1.00 40.36           C  
ATOM   2677  CZ  PHE A 348       9.635  -5.187  31.265  1.00 39.92           C  
""",
    {
        Feature('AROM', [10.44, -1.877, 29.88]),
        Feature('AROM', [11.06, -5.208, 31.38]),
        Feature('AROM', [11.64, -4.319, 32.37]),
        Feature('AROM', [8.773, -4.363, 32.04]),
        Feature('AROM', [9.374, -3.485, 33.03]),
        Feature('AROM', [9.607, -5.225, 31.22]),
        Feature('AROM', [9.98, -6.801, 34.54]),
        Feature('HACC', [11.3, -1.116, 36.78]),
        Feature('HDON', [13.45, -5.081, 36.05]),
        Feature('LIPO', [10.16, -4.919, 32.76]),
        Feature('LIPO', [10.27, -3.76, 31.66]),
        Feature('LIPO', [12.16, -6.24, 30.4]),
        Feature('LIPO', [13.43, -4.316, 32.55]),
        Feature('LIPO', [6.982, -4.366, 31.86]),
        Feature('LIPO', [8.298, -2.441, 34.02]),
        Feature('LIPO', [8.871, -6.318, 29.99]),
    },
    id='PHE'
)

PRO = pytest.param(
    """ATOM   2692  N   PRO A 351       8.691  -9.736  32.054  1.00 54.28           N  
ATOM   2693  CA  PRO A 351       7.695 -10.810  32.164  1.00 55.40           C  
ATOM   2694  C   PRO A 351       8.256 -12.172  31.720  1.00 56.38           C  
ATOM   2695  O   PRO A 351       9.265 -12.216  30.997  1.00 56.66           O  
ATOM   2696  CB  PRO A 351       6.581 -10.359  31.217  1.00 55.19           C  
ATOM   2697  CG  PRO A 351       7.292  -9.516  30.195  1.00 55.29           C  
ATOM   2698  CD  PRO A 351       8.392  -8.815  30.940  1.00 54.62           C  
""",
    {
        Feature('HDON', [9.879, -12.31, 30.5]),
        Feature('LIPO', [5.502, -9.168, 32.26]),
        Feature('LIPO', [5.633, -11.82, 30.5]),
        Feature('LIPO', [6.149, -8.343, 29.25]),
        Feature('LIPO', [7.637, -7.251, 31.73]),
        Feature('LIPO', [8.172, -10.74, 28.99]),
        Feature('LIPO', [9.956, -8.49, 29.94]),
    },
    id='PRO'
)

SER = pytest.param(
    """ATOM   2661  N   SER A 347      13.282  -0.315  38.196  1.00 41.41           N  
ATOM   2662  CA  SER A 347      12.393  -1.290  38.846  1.00 41.94           C  
ATOM   2663  C   SER A 347      11.587  -2.191  37.914  1.00 42.47           C  
ATOM   2664  O   SER A 347      10.670  -2.885  38.369  1.00 42.91           O  
ATOM   2665  CB  SER A 347      11.431  -0.576  39.783  1.00 42.38           C  
ATOM   2666  OG  SER A 347      10.474   0.173  39.047  1.00 44.62           O  
""",
    {
        Feature('HACC', [13.81, -0.6285, 37.3]),
        Feature('HDON', [11.44, 0.463, 39.8]),
        Feature('HDON', [9.842, -3.037, 38.06]),
    },
    id='SER'
)

THR = pytest.param(
    """ATOM   2610  N   THR A 341      10.795   8.648  36.299  1.00 30.96           N  
ATOM   2611  CA  THR A 341      10.853   7.768  35.112  1.00 29.85           C  
ATOM   2612  C   THR A 341      12.081   6.831  35.213  1.00 31.19           C  
ATOM   2613  O   THR A 341      11.986   5.608  34.919  1.00 30.11           O  
ATOM   2614  CB  THR A 341      10.833   8.581  33.788  1.00 29.22           C  
ATOM   2615  OG1 THR A 341       9.600   9.344  33.701  1.00 24.89           O  
ATOM   2616  CG2 THR A 341      10.955   7.665  32.565  1.00 27.62           C  
""",
    {
        Feature('HACC', [9.604, 9.076, 36.25]),
        Feature('HDON', [12.02, 4.849, 34.74]),
        Feature('HDON', [9.799, 8.664, 32.97]),
        Feature('LIPO', [10.9, 7.13, 31.95]),
        Feature('LIPO', [11.73, 8.508, 31.1]),
        Feature('LIPO', [11.8, 6.021, 32.96]),
        Feature('LIPO', [9.047, 7.241, 32.22]),
    },
    id='THR'
)

TRP = pytest.param(
    """ATOM   2573  N   TRP A 337      10.077  14.558  37.421  1.00 27.64           N  
ATOM   2574  CA  TRP A 337       9.446  13.667  36.455  1.00 27.45           C  
ATOM   2575  C   TRP A 337      10.464  12.735  35.847  1.00 28.20           C  
ATOM   2576  O   TRP A 337      10.206  11.540  35.687  1.00 29.47           O  
ATOM   2577  CB  TRP A 337       8.757  14.463  35.352  1.00 26.13           C  
ATOM   2578  CG  TRP A 337       7.490  15.137  35.810  1.00 25.39           C  
ATOM   2579  CD1 TRP A 337       7.161  16.449  35.666  1.00 24.76           C  
ATOM   2580  CD2 TRP A 337       6.381  14.516  36.485  1.00 25.73           C  
ATOM   2581  NE1 TRP A 337       5.919  16.692  36.219  1.00 24.24           N  
ATOM   2582  CE2 TRP A 337       5.423  15.521  36.728  1.00 25.45           C  
ATOM   2583  CE3 TRP A 337       6.106  13.196  36.908  1.00 27.13           C  
ATOM   2584  CZ2 TRP A 337       4.198  15.255  37.355  1.00 25.01           C  
ATOM   2585  CZ3 TRP A 337       4.912  12.946  37.559  1.00 25.10           C  
ATOM   2586  CH2 TRP A 337       3.968  13.973  37.767  1.00 23.56           C  
    """,
    {
        Feature('AROM', [10.07, 14.58, 37.37]),
        Feature('AROM', [10.12, 14.59, 37.32]),
        Feature('AROM', [10.52, 12.75, 35.86]),
        Feature('AROM', [5.919, 16.73, 36.25]),
        Feature('AROM', [6.086, 13.01, 33.7]),
        Feature('AROM', [6.653, 13.94, 33.17]),
        Feature('AROM', [7.208, 16.5, 35.66]),
        Feature('AROM', [8.703, 14.33, 35.26]),
        Feature('AROM', [8.786, 14.48, 35.32]),
        Feature('AROM', [9.052, 14.18, 39.7]),
        Feature('AROM', [9.415, 13.64, 36.53]),
        Feature('AROM', [9.62, 15.12, 39.17]),
        Feature('HACC', [5.355, 17.79, 36.33]),
        Feature('HACC', [9.508, 14.98, 38.28]),
        Feature('HDON', [10.06, 10.77, 35.61]),
        Feature('LIPO', [11.48, 15.49, 36.72]),
        Feature('LIPO', [12.14, 13.37, 35.37]),
        Feature('LIPO', [5.073, 18.32, 36.36]),
        Feature('LIPO', [7.788, 14.39, 35.46]),
        Feature('LIPO', [8.223, 17.78, 34.91]),
        Feature('LIPO', [8.251, 13.12, 34.01]),
        Feature('LIPO', [8.36, 12.63, 37.58]),
        Feature('LIPO', [8.486, 14.67, 36.88]),
        Feature('LIPO', [9.202, 15.17, 38.75]),
        Feature('LIPO', [9.896, 15.76, 34.72]),
    },
    id='TRP'
)

TYR = pytest.param(
    """ATOM     37  N   TYR A   9      24.803  -1.481  21.129  1.00 55.74           N  
ATOM     38  CA  TYR A   9      25.410  -1.127  19.863  1.00 56.50           C  
ATOM     39  C   TYR A   9      24.669   0.056  19.256  1.00 57.38           C  
ATOM     40  O   TYR A   9      23.440   0.124  19.319  1.00 57.50           O  
ATOM     41  CB  TYR A   9      25.486  -2.340  18.917  1.00 56.03           C  
ATOM     42  CG  TYR A   9      24.164  -2.873  18.402  1.00 55.82           C  
ATOM     43  CD1 TYR A   9      23.739  -2.581  17.103  1.00 56.44           C  
ATOM     44  CD2 TYR A   9      23.350  -3.682  19.194  1.00 54.64           C  
ATOM     45  CE1 TYR A   9      22.533  -3.067  16.613  1.00 55.11           C  
ATOM     46  CE2 TYR A   9      22.146  -4.174  18.712  1.00 54.00           C  
ATOM     47  CZ  TYR A   9      21.745  -3.861  17.422  1.00 55.25           C  
ATOM     48  OH  TYR A   9      20.553  -4.329  16.923  1.00 54.78           O  
    """,
    {
        Feature('AROM', [20.56, -4.357, 16.91]),
        Feature('AROM', [22.47, 1.255, 20.4]),
        Feature('AROM', [24.71, 0.1081, 19.23]),
        Feature('AROM', [24.81, -1.478, 21.11]),
        Feature('AROM', [24.83, -1.493, 21.11]),
        Feature('AROM', [25.49, -1.093, 19.88]),
        Feature('AROM', [25.52, -2.406, 18.97]),
        Feature('AROM', [25.54, -2.317, 18.84]),
        Feature('AROM', [25.65, -4.288, 18.07]),
        Feature('HACC', [23.89, -0.8354, 21.45]),
        Feature('HDON', [21.3, -4.911, 16.64]),
        Feature('HDON', [22.65, 0.179, 19.34]),
        Feature('LIPO', [20.64, -5.854, 15.91]),
        Feature('LIPO', [23.05, -1.592, 18.96]),
        Feature('LIPO', [23.41, -0.5067, 21.61]),
        Feature('LIPO', [23.8, -2.897, 18.41]),
        Feature('LIPO', [25.63, 1.42, 18.42]),
        Feature('LIPO', [25.95, -2.239, 22.28]),
        Feature('LIPO', [26.24, -3.711, 19.98]),
        Feature('LIPO', [26.55, -1.849, 17.43]),
        Feature('LIPO', [27.21, -0.598, 20.03])
    },
    id='TYR'
)

# 3HEG has a TYR on A 35 which is truncated
TYR_TRUNCATED = pytest.param(
    """ATOM    232  N   TYR A  35       3.890   5.530  15.364  1.00 56.36           N  
ATOM    233  CA  TYR A  35       4.070   6.446  16.490  1.00 55.91           C  
ATOM    234  C   TYR A  35       5.048   7.601  16.179  1.00 55.30           C  
ATOM    235  O   TYR A  35       4.746   8.521  15.397  1.00 55.26           O  
ATOM    236  CB  TYR A  35       4.508   5.658  17.762  1.00 55.74           C  
    """,
    {
        Feature('HACC', [3.94, 4.382, 15.59]),
        Feature('HDON', [4.596, 9.111, 14.88]),
    },
    id='TYR_TRUNCATED'
)

VAL = pytest.param(
    """ATOM    225  N   VAL A  30       9.593    1.42  11.006  1.00 52.99           N
ATOM    226  CA  VAL A  30       8.249   0.868  10.857  1.00 53.35           C
ATOM    227  C   VAL A  30       7.187   1.933  11.164  1.00 54.25           C
ATOM    228  O   VAL A  30       6.129   2.018  10.498  1.00 54.85           O
ATOM    229  CB  VAL A  30       8.097  -0.423  11.743  1.00  53.8           C
ATOM    230  CG1 VAL A  30       7.479  -0.127  13.096  1.00 52.84           C
ATOM    231  CG2 VAL A  30       7.324  -1.515  11.024  1.00 53.61           C
    """,
    {
        Feature('HACC', [10.18, 1.059, 11.95]),
        Feature('HDON', [5.626, 2.703, 11.3]),
        Feature('LIPO', [6.95, 1.571, 13.53]),
        Feature('LIPO', [8.149, 0.07703, 14.65]),
        Feature('LIPO', [6.158, -0.4736, 13.86]),
        Feature('LIPO', [6.783, 1.71, 13.02]),
        Feature('LIPO', [6.969, -2.094, 10.67]),
        Feature('LIPO', [5.8, -1.179, 13.25]),
        Feature('LIPO', [9.937, -0.9608, 12.11]),
        Feature('LIPO', [6.978, 0.05671, 13.71]),
        Feature('LIPO', [8.443, -0.3411, 14.58]),
    },
    id='VAL'
)


def feat2point(f: Feature):
    return f.kind, f.position[0], f.position[1], f.position[2]


def feat2points(features):
    return [feat2point(f) for f in features]


@pytest.mark.parametrize("block,expected", [
    # ALA,
    # ARG,
    # ASN,
    # ASP,
    # CYS,
    # HIS,
    # GLU,
    # GLN,
    # GLY,
    ILE,
    # LEU,
    # LYS,
    # MET,
    # PHE,
    # PRO,
    # SER,
    # THR,
    # TRP,
    # TYR,
    # TYR_TRUNCATED,
    # VAL,
])
def test_from_site(block, expected, ligand_3heg_bax: Ligand):
    site = prep_site(block, ligand_3heg_bax.molecule)

    features = from_site(site)
    # print(site.ligand().model().to_file_string('pdb'))
    print(features)
    # dump4molviewer(features, site)

    assert_features(expected, features)


def dump4molviewer(features, site):
    label = site.residue().name()
    data = [{
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
    }]
    with open('ALA.json', 'w') as f:
        json.dump(data, f)
