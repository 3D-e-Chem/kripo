from typing import Set

from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.structures.chains import Site
import pytest
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem.rdmolfiles import MolFromPDBBlock, MolToPDBBlock
from rdkit.Chem.rdmolops import AddHs

from kripo.pharmacophore import Feature, from_site


def prep_site(block):
    mol = MolFromPDBBlock(block)
    # Hydronate molecule
    m2 = AddHs(mol)
    coordMap = {}
    atomMap = []
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        coordMap[i] = conf.GetAtomPosition(i)
        atomMap.append((i, i,))

    # params = ETKDG()
    # params.coordMap = coordMap
    # params.randomSeed = 42
    # EmbedMolecule(m2, params)
    EmbedMolecule(m2, coordMap=coordMap, randomSeed=42, useBasicKnowledge=True)
    AlignMol(m2, mol, atomMap=atomMap)

    protonated_block = MolToPDBBlock(m2)
    ligand = """HETATM 2706  C1  BAX A   1      -1.978   2.259  25.948  1.00 25.24           C
"""
    out_block = ligand + protonated_block
    model = pdb_dict_to_pdb(pdb_string_to_pdb_dict(out_block)).model()
    return Site(model, ligand=model.molecule(name='BAX'))


def assert_features(expected: Set[Feature], result: Set[Feature]):
    sorted_expected = sorted(list(expected), key=lambda d: repr(d))
    sorted_result = sorted(list(result), key=lambda d: repr(d))
    assert sorted_expected == sorted_result


ALA = pytest.param(
    """ATOM   2432  N   ALA A 320      -7.954  -5.611  36.611  1.00 37.45           N  
ATOM   2433  CA  ALA A 320      -6.613  -5.159  36.864  1.00 37.02           C  
ATOM   2434  C   ALA A 320      -6.315  -4.998  38.349  1.00 37.13           C  
ATOM   2435  O   ALA A 320      -7.194  -4.648  39.136  1.00 35.74           O  
ATOM   2436  CB  ALA A 320      -6.373  -3.834  36.140  1.00 37.00           C  
""",
    {
        Feature('HACC', [-8.085, -6.55, 35.93]),
        Feature('HDON', [-7.758, -4.466, 39.67]),
        Feature('LIPO', [-4.498, -3.721, 35.86]),
        Feature('LIPO', [-6.233, -3.152, 35.75]),
        Feature('LIPO', [-6.983, -2.367, 37.18]),
        Feature('LIPO', [-7.298, -3.982, 34.47]),
    },
    id='ALA'
)

ARG = pytest.param(
    """ATOM   2513  N   ARG A 330      -0.182  16.030  40.538  1.00 35.62           N  
ATOM   2514  CA  ARG A 330       1.082  16.624  40.900  1.00 34.90           C  
ATOM   2515  C   ARG A 330       1.609  17.501  39.794  1.00 34.79           C  
ATOM   2516  O   ARG A 330       1.643  17.107  38.629  1.00 35.15           O  
ATOM   2517  CB  ARG A 330       2.087  15.519  41.179  1.00 34.43           C  
ATOM   2518  CG  ARG A 330       2.457  15.373  42.611  1.00 36.84           C  
ATOM   2519  CD  ARG A 330       1.431  14.636  43.422  1.00 38.13           C  
ATOM   2520  NE  ARG A 330       1.783  13.230  43.454  1.00 41.26           N  
ATOM   2521  CZ  ARG A 330       1.226  12.328  44.248  1.00 42.35           C  
ATOM   2522  NH1 ARG A 330       0.293  12.671  45.128  1.00 42.54           N  
ATOM   2523  NH2 ARG A 330       1.616  11.072  44.163  1.00 40.54           N  
""",
    {
        Feature('HACC', [-0.7258, 16.18, 41.84]),
        Feature('HACC', [-0.8002, 16.4, 39.85]),
        Feature('HACC', [-0.827, 11.63, 44.49]),
        Feature('HACC', [0.7588, 10.73, 43.45]),
        Feature('HACC', [1.067, 14.8, 44.09]),
        Feature('HACC', [2.79, 13.27, 45.98]),
        Feature('HDON', [2.258, 18.93, 39.81]),
        Feature('LIPO', [-0.2348, 14.15, 41.04]),
        Feature('LIPO', [2.177, 13.02, 39.71]),
        Feature('LIPO', [3.304, 16.01, 42.59]),
        Feature('LIPO', [3.468, 15.48, 39.74]),
        Feature('POSC', [-0.2318, 16.56, 40.62]),
        Feature('POSC', [0.2163, 11.48, 44.09]),
        Feature('POSC', [1.498, 13.85, 43.54]),
        Feature('POSC', [1.858, 12.64, 46.06]),
    },
    id='ARG'
)

ASN = pytest.param(
    """ATOM   2102  N   ASN A 278     -27.300  -9.279  12.569  1.00 43.44           N  
ATOM   2103  CA  ASN A 278     -27.866 -10.088  13.617  1.00 42.99           C  
ATOM   2104  C   ASN A 278     -28.685  -9.180  14.542  1.00 43.56           C  
ATOM   2105  O   ASN A 278     -28.167  -8.162  15.003  1.00 43.72           O  
ATOM   2106  CB  ASN A 278     -26.746 -10.798  14.379  1.00 42.39           C  
ATOM   2107  CG  ASN A 278     -27.250 -11.667  15.501  1.00 41.27           C  
ATOM   2108  OD1 ASN A 278     -28.381 -11.541  15.954  1.00 41.81           O  
ATOM   2109  ND2 ASN A 278     -26.402 -12.558  15.970  1.00 42.33           N  
""",
    {
        Feature('HACC', [-25.24, -12.66, 15.84]),
        Feature('HACC', [-26.23, -8.799, 12.74]),
        Feature('HACC', [-26.95, -13.76, 16.34]),
        Feature('HDON', [-28.98, -11.27, 16.45]),
        Feature('HDON', [-29.67, -8.153, 14.76]),
    },
    id='ASN'
)

ASP = pytest.param(
    """ATOM   2629  N   ASP A 343      13.634   5.771  37.980  1.00 32.80           N  
ATOM   2630  CA  ASP A 343      13.219   4.716  38.898  1.00 34.48           C  
ATOM   2631  C   ASP A 343      12.388   3.625  38.203  1.00 34.84           C  
ATOM   2632  O   ASP A 343      12.693   2.426  38.327  1.00 36.08           O  
ATOM   2633  CB  ASP A 343      12.395   5.291  40.036  1.00 34.50           C  
ATOM   2634  CG  ASP A 343      13.158   6.301  40.844  1.00 37.06           C  
ATOM   2635  OD1 ASP A 343      12.493   7.082  41.578  1.00 40.45           O  
ATOM   2636  OD2 ASP A 343      14.409   6.315  40.726  1.00 36.59           O  
""",
    {
        Feature('HACC', [14.2, 5.575, 37.34]),
        Feature('HDON', [12.0, 1.832, 38.22]),
        Feature('HDON', [13.17, 7.841, 39.33]),
        Feature('HDON', [13.29, 6.682, 42.79]),
        Feature('NEGC', [12.25, 8.483, 39.01]),
        Feature('NEGC', [12.41, 6.983, 43.48]),
        Feature('NEGC', [13.35, 7.644, 41.31]),
        Feature('NEGC', [14.19, 8.045, 38.8]),
        Feature('NEGC', [14.35, 6.545, 43.27]),
    },
    id='ASP'
)

CYS = pytest.param(
    """ATOM   1578  N   CYS A 211     -19.537   2.552  17.945  1.00 29.42           N  
ATOM   1579  CA  CYS A 211     -19.222   3.112  16.636  1.00 29.59           C  
ATOM   1580  C   CYS A 211     -18.394   2.164  15.742  1.00 30.27           C  
ATOM   1581  O   CYS A 211     -18.610   2.103  14.524  1.00 30.26           O  
ATOM   1582  CB  CYS A 211     -18.470   4.439  16.799  1.00 29.76           C  
ATOM   1583  SG  CYS A 211     -19.482   5.843  17.432  1.00 29.04           S  
""",
    {
        Feature('HACC', [-18.64, 2.087, 18.53]),
        Feature('HDON', [-18.69, 2.055, 13.75]),
        Feature('LIPO', [-21.55, 5.366, 17.76]),
    },
    id='CYS'
)

HIS = pytest.param(
    """ATOM   2368  N   HIS A 312     -18.230 -15.709  25.729  1.00 41.27           N  
ATOM   2369  CA  HIS A 312     -18.266 -15.219  27.094  1.00 41.07           C  
ATOM   2370  C   HIS A 312     -17.913 -16.336  28.043  1.00 41.80           C  
ATOM   2371  O   HIS A 312     -18.423 -17.441  27.919  1.00 41.92           O  
ATOM   2372  CB  HIS A 312     -19.660 -14.664  27.389  1.00 40.71           C  
ATOM   2373  CG  HIS A 312     -19.882 -14.217  28.808  1.00 41.18           C  
ATOM   2374  ND1 HIS A 312     -19.023 -13.373  29.475  1.00 39.83           N  
ATOM   2375  CD2 HIS A 312     -20.919 -14.438  29.654  1.00 40.93           C  
ATOM   2376  CE1 HIS A 312     -19.500 -13.127  30.683  1.00 41.50           C  
ATOM   2377  NE2 HIS A 312     -20.648 -13.765  30.816  1.00 43.14           N  
""",
    {
        Feature('AROM', [-16.87, -18.61, 28.44]),
        Feature('AROM', [-20.13, -13.15, 26.03]),
        Feature('HACC', [-17.1, -15.86, 25.3]),
        Feature('HACC', [-21.37, -13.75, 31.83]),
        Feature('HDON', [-18.74, -18.17, 27.85]),
    },
    id='HIS'
)

GLU = pytest.param(
    """ATOM   2637  N   GLU A 344      11.332   4.036  37.511  1.00 33.61           N  
ATOM   2638  CA  GLU A 344      10.446   3.091  36.858  1.00 34.11           C  
ATOM   2639  C   GLU A 344      11.193   2.204  35.849  1.00 34.20           C  
ATOM   2640  O   GLU A 344      10.876   1.033  35.686  1.00 34.31           O  
ATOM   2641  CB  GLU A 344       9.251   3.810  36.226  1.00 33.38           C  
ATOM   2642  CG  GLU A 344       8.284   4.395  37.233  1.00 34.10           C  
ATOM   2643  CD  GLU A 344       7.540   3.341  38.086  1.00 35.64           C  
ATOM   2644  OE1 GLU A 344       7.358   2.184  37.644  1.00 34.95           O  
ATOM   2645  OE2 GLU A 344       7.078   3.704  39.188  1.00 35.93           O  
""",
    {
        Feature('HACC', [11.99, 4.744, 36.83]),
        Feature('HDON', [10.66, 0.2883, 35.5]),
        Feature('HDON', [6.792, 3.91, 39.91]),
        Feature('HDON', [7.252, 1.434, 37.38]),
        Feature('NEGC', [5.718, 4.351, 40.01]),
        Feature('NEGC', [6.311, 1.153, 36.75]),
        Feature('NEGC', [6.844, 2.398, 38.82]),
        Feature('NEGC', [7.436, 3.792, 40.87]),
        Feature('NEGC', [8.028, 0.5946, 37.61]),
    },
    id='GLU'
)

GLN = pytest.param(
    """ATOM     60  N   GLN A  11      23.638   2.737  16.073  1.00 61.44           N  
ATOM     61  CA  GLN A  11      23.420   2.838  14.623  1.00 62.58           C  
ATOM     62  C   GLN A  11      22.559   4.046  14.201  1.00 63.43           C  
ATOM     63  O   GLN A  11      21.818   4.615  15.004  1.00 63.68           O  
ATOM     64  CB  GLN A  11      22.868   1.532  14.031  1.00 62.37           C  
ATOM     65  CG  GLN A  11      21.560   1.065  14.625  1.00 62.67           C  
ATOM     66  CD  GLN A  11      20.972  -0.108  13.865  1.00 62.48           C  
ATOM     67  OE1 GLN A  11      20.666  -0.003  12.676  1.00 62.75           O  
ATOM     68  NE2 GLN A  11      20.810  -1.233  14.549  1.00 61.50           N  
""",
    {
        Feature('HACC', [20.64, -1.309, 15.8]),
        Feature('HACC', [20.78, -2.376, 14.0]),
        Feature('HACC', [23.06, 3.421, 16.78]),
        Feature('HDON', [20.48, 0.06548, 11.9]),
        Feature('HDON', [21.31, 5.013, 15.47]),
    },
    id='GLN'
)

GLY = pytest.param(
    """ATOM   2093  N   GLY A 276     -30.316  -6.295   7.777  1.00 49.02           N  
ATOM   2094  CA  GLY A 276     -30.341  -7.610   8.362  1.00 47.92           C  
ATOM   2095  C   GLY A 276     -28.954  -7.961   8.885  1.00 47.11           C  
ATOM   2096  O   GLY A 276     -28.325  -8.879   8.375  1.00 47.06           O  
""",
    {
        Feature('HACC', [-29.16, -5.804, 7.884]),
        Feature('HDON', [-27.86, -9.447, 8.052]),
    },
    id='GLY'
)

ILE = pytest.param(
    """ATOM   2653  N   ILE A 346      14.387   1.339  36.234  1.00 39.47           N  
ATOM   2654  CA  ILE A 346      15.274   0.501  37.045  1.00 40.11           C  
ATOM   2655  C   ILE A 346      14.485  -0.631  37.720  1.00 40.36           C  
ATOM   2656  O   ILE A 346      14.947  -1.772  37.758  1.00 39.90           O  
ATOM   2657  CB  ILE A 346      16.118   1.379  38.045  1.00 40.58           C  
ATOM   2658  CG1 ILE A 346      17.424   1.872  37.408  1.00 41.02           C  
ATOM   2659  CG2 ILE A 346      16.478   0.625  39.358  1.00 40.14           C  
ATOM   2660  CD1 ILE A 346      17.601   1.631  35.886  1.00 43.16           C  
""",
    {
        Feature('HACC', [13.93, 0.8756, 35.3]),
        Feature('HDON', [14.28, -1.752, 38.93]),
        Feature('LIPO', [14.89, 2.579, 38.75]),
        Feature('LIPO', [15.52, 0.8621, 40.78]),
        Feature('LIPO', [16.13, 2.268, 34.89]),
        Feature('LIPO', [16.76, -1.361, 38.85]),
        Feature('LIPO', [16.86, 0.07399, 39.91]),
        Feature('LIPO', [17.22, 3.866, 37.73]),
        Feature('LIPO', [17.73, 1.609, 35.16]),
        Feature('LIPO', [18.23, -0.03753, 35.67]),
        Feature('LIPO', [18.38, 1.055, 39.71]),
        Feature('LIPO', [18.96, 1.423, 38.29]),
        Feature('LIPO', [19.04, 2.887, 35.65]),
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


@pytest.mark.parametrize("block,expected", [
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
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
])
def test_from_site(block, expected):
    site = prep_site(block)

    features = from_site(site)

    assert_features(expected, features)
