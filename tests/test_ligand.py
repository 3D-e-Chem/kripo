from kripo.ligand import Ligand, remove_nonpdb_bonds
from rdkit.Chem.rdmolfiles import MolFromPDBBlock, MolToSmiles
from rdkit.Chem.rdmolops import SanitizeMol
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict

def test_name(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.name() == 'BAX'


def test_pdb_block(ligand_3heg_bax: Ligand):
    block = ligand_3heg_bax.pdb_block()
    assert 'HETATM' in block
    assert 'BAX' in block
    assert 'CL11' in block
    assert 'CONECT' in block


def test_fragments(ligand_3heg_bax: Ligand):
    fragments = ligand_3heg_bax.fragments()

    assert len(fragments) == 28


def test_id(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.id() == 'A1'


def test_chain(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.chain() == 'A'


def test_seq_nr(ligand_3heg_bax: Ligand):
    assert ligand_3heg_bax.seq_nr() == 1


def test_remove_nonpdb_bonds__norendundancy_returnself():
    block = """HETATM 2835  O2  SB2 A 800      -0.940  21.570  27.801  1.00 75.60           O  
HETATM 2836  C1  SB2 A 800       1.767  22.346  28.287  1.00 75.48           C  
HETATM 2837  S1  SB2 A 800       0.588  21.264  27.403  1.00 75.07           S  
HETATM 2838  CA1 SB2 A 800       0.463  17.215  27.773  1.00 66.70           C  
HETATM 2839  CA2 SB2 A 800       0.075  18.563  27.590  1.00 69.49           C  
HETATM 2840  CA3 SB2 A 800       1.069  19.579  27.613  1.00 71.50           C  
HETATM 2841  CA4 SB2 A 800       2.444  19.278  27.802  1.00 68.18           C  
HETATM 2842  CA5 SB2 A 800       2.826  17.933  27.981  1.00 66.56           C  
HETATM 2843  CA6 SB2 A 800       1.852  16.853  27.972  1.00 63.47           C  
HETATM 2844  NB1 SB2 A 800       0.469   9.354  28.905  1.00 51.18           N  
HETATM 2845  CB2 SB2 A 800       0.033  10.130  27.884  1.00 51.09           C  
HETATM 2846  CB3 SB2 A 800       0.560  11.449  27.685  1.00 51.52           C  
HETATM 2847  CB4 SB2 A 800       1.572  11.953  28.598  1.00 54.09           C  
HETATM 2848  CB5 SB2 A 800       1.989  11.057  29.683  1.00 52.94           C  
HETATM 2849  CB6 SB2 A 800       1.408   9.765  29.789  1.00 51.35           C  
HETATM 2850  NC1 SB2 A 800       1.347  14.432  28.206  1.00 55.23           N  
HETATM 2851  CC2 SB2 A 800       2.255  15.490  28.183  1.00 57.53           C  
HETATM 2852  NC3 SB2 A 800       3.576  15.169  28.381  1.00 52.86           N  
HETATM 2853  CC4 SB2 A 800       3.553  13.784  28.563  1.00 52.74           C  
HETATM 2854  CC5 SB2 A 800       2.150  13.303  28.479  1.00 53.92           C  
HETATM 2855  CD1 SB2 A 800       5.338  13.337  30.169  1.00 50.85           C  
HETATM 2856  CD2 SB2 A 800       6.535  12.707  30.543  1.00 52.79           C  
HETATM 2857  CD3 SB2 A 800       7.189  11.834  29.638  1.00 51.62           C  
HETATM 2858  CD4 SB2 A 800       6.636  11.584  28.339  1.00 49.80           C  
HETATM 2859  CD5 SB2 A 800       5.427  12.219  27.960  1.00 49.08           C  
HETATM 2860  CD6 SB2 A 800       4.755  13.110  28.880  1.00 51.23           C  
HETATM 2861  FD3 SB2 A 800       8.355  11.233  30.035  1.00 53.01           F  
CONECT 2835 2837                                                                
CONECT 2836 2837                                                                
CONECT 2837 2835 2836 2840                                                      
CONECT 2838 2839 2843                                                           
CONECT 2839 2838 2840                                                           
CONECT 2840 2837 2839 2841                                                      
CONECT 2841 2840 2842                                                           
CONECT 2842 2841 2843                                                           
CONECT 2843 2838 2842 2851                                                      
CONECT 2844 2845 2849       
CONECT 2845 2844 2846                                                           
CONECT 2846 2845 2847                                                           
CONECT 2847 2846 2848 2854                                                      
CONECT 2848 2847 2849                                                           
CONECT 2849 2844 2848                                                           
CONECT 2850 2851 2854                                                           
CONECT 2851 2843 2850 2852                                                      
CONECT 2852 2851 2853                                                           
CONECT 2853 2852 2854 2860                                                      
CONECT 2854 2847 2850 2853                                                      
CONECT 2855 2856 2860                                                           
CONECT 2856 2855 2857                                                           
CONECT 2857 2856 2858 2861                                                      
CONECT 2858 2857 2859                                                           
CONECT 2859 2858 2860                                                           
CONECT 2860 2853 2855 2859                                                      
CONECT 2861 2857     
"""

    rmol = MolFromPDBBlock(block)
    amol = pdb_dict_to_pdb(pdb_string_to_pdb_dict(block)).model().molecule()

    result = remove_nonpdb_bonds(rmol, amol)

    assert MolToSmiles(rmol) == MolToSmiles(result)


def test_remove_nonpdb_bonds__withredundancy_sanitizable():
    block = """HETATM 2835  O2  SB2 A 800      -0.940  21.570  27.801  1.00  75.6           O  
HETATM 2836  C1  SB2 A 800       1.767  22.346  28.287  1.00 75.48           C  
HETATM 2837  S1  SB2 A 800       0.588  21.264  27.403  1.00 75.07           S  
HETATM 2838  CA1 SB2 A 800       0.463  17.215  27.773  1.00  66.7           C  
HETATM 2839  CA2 SB2 A 800       0.075  18.563  27.590  1.00 69.49           C  
HETATM 2840  CA3 SB2 A 800       1.069  19.579  27.613  1.00  71.5           C  
HETATM 2841  CA4 SB2 A 800       2.444  19.278  27.802  1.00 68.18           C  
HETATM 2842  CA5 SB2 A 800       2.826  17.933  27.981  1.00 66.56           C  
HETATM 2843  CA6 SB2 A 800       1.852  16.853  27.972  1.00 63.47           C  
HETATM 2844  NB1 SB2 A 800       0.469   9.354  28.905  1.00 51.18           N  
HETATM 2845  CB2 SB2 A 800       0.033  10.130  27.884  1.00 51.09           C  
HETATM 2846  CB3 SB2 A 800       0.560  11.449  27.685  1.00 51.52           C  
HETATM 2847  CB4 SB2 A 800       1.572  11.953  28.598  1.00 54.09           C  
HETATM 2848  CB5 SB2 A 800       1.989  11.057  29.683  1.00 52.94           C  
HETATM 2849  CB6 SB2 A 800       1.408   9.765  29.789  1.00 51.35           C  
HETATM 2850  NC1 SB2 A 800       1.347  14.432  28.206  1.00 55.23           N  
HETATM 2851  CC2 SB2 A 800       2.255  15.490  28.183  1.00 57.53           C  
HETATM 2852  NC3 SB2 A 800       3.576  15.169  28.381  1.00 52.86           N  
HETATM 2853  CC4 SB2 A 800       3.553  13.784  28.563  1.00 52.74           C  
HETATM 2854  CC5 SB2 A 800       2.150  13.303  28.479  1.00 53.92           C  
HETATM 2855  CD1 SB2 A 800       5.338  13.337  30.169  1.00 50.85           C  
HETATM 2856  CD2 SB2 A 800       6.535  12.707  30.543  1.00 52.79           C  
HETATM 2857  CD3 SB2 A 800       7.189  11.834  29.638  1.00 51.62           C  
HETATM 2858  CD4 SB2 A 800       6.636  11.584  28.339  1.00  49.8           C  
HETATM 2859  CD5 SB2 A 800       5.427  12.219  27.960  1.00 49.08           C  
HETATM 2860  CD6 SB2 A 800       4.755  13.110  28.880  1.00 51.23           C  
HETATM 2861  FD3 SB2 A 800       8.355  11.233  30.035  1.00 53.01           F  
HETATM 2862  H   SB2 A 800       4.870  13.974  30.833  1.00                 H  
HETATM 2863  H   SB2 A 800       3.151  20.030  27.808  1.00                 H  
HETATM 2864  H   SB2 A 800       3.823  17.709  28.122  1.00                 H  
HETATM 2865  H   SB2 A 800       2.232  12.273  28.758  1.00                 H  
HETATM 2866  H   SB2 A 800      -0.688   9.767  27.241  1.00                 H  
HETATM 2867  H   SB2 A 800       6.939  12.879  31.477  1.00                 H  
HETATM 2868  H   SB2 A 800       5.025  12.045  27.025  1.00                 H  
HETATM 2869  H   SB2 A 800       7.115  10.946  27.684  1.00                 H  
HETATM 2870  H   SB2 A 800       0.225  12.033  26.903  1.00                 H  
HETATM 2871  H   SB2 A 800       1.705   9.133  30.549  1.00                 H  
HETATM 2872  H   SB2 A 800       2.750  22.203  27.889  1.00                 H  
HETATM 2873  H   SB2 A 800       1.477  23.368  28.161  1.00                 H  
HETATM 2874  H   SB2 A 800       1.764  22.101  29.328  1.00                 H  
HETATM 2875  H   SB2 A 800       2.700  11.364  30.365  1.00                 H  
HETATM 2876  H   SB2 A 800      -0.917  18.806  27.442  1.00                 H  
HETATM 2877  H   SB2 A 800      -0.257  16.475  27.765  1.00                 H  
CONECT 2835 2837                                                                
CONECT 2836 2837 2872 2873 2874                                                 
CONECT 2837 2835 2836 2840                                                      
CONECT 2838 2839 2843 2877                                                      
CONECT 2839 2838 2840 2876                                                      
CONECT 2840 2837 2839 2841                                                      
CONECT 2841 2840 2842 2863                                                      
CONECT 2842 2841 2843 2864                                                      
CONECT 2843 2838 2842 2851                                                      
CONECT 2844 2845 2849                                                           
CONECT 2845 2844 2846 2866                                                      
CONECT 2846 2845 2847 2870                                                      
CONECT 2847 2846 2848 2854                                                      
CONECT 2848 2847 2849 2875                                                      
CONECT 2849 2844 2848 2871                                                      
CONECT 2850 2851 2854                                                           
CONECT 2851 2843 2850 2852                                                      
CONECT 2852 2851 2853                                                           
CONECT 2853 2852 2854 2860                                                      
CONECT 2854 2847 2850 2853 2865                                                 
CONECT 2855 2856 2860 2862                                                      
CONECT 2856 2855 2857 2867                                                      
CONECT 2857 2856 2858 2861                                                      
CONECT 2858 2857 2859 2869                                                      
CONECT 2859 2858 2860 2868                                                      
CONECT 2860 2853 2855 2859                                                      
CONECT 2861 2857                                                                
CONECT 2862 2855                                                                
CONECT 2863 2841                                                                
CONECT 2864 2842                                                                
CONECT 2865 2854                                                                
CONECT 2866 2845                                                                
CONECT 2867 2856                                                                
CONECT 2868 2859                                                                
CONECT 2869 2858                                                                
CONECT 2870 2846                                                                
CONECT 2871 2849                                                                
CONECT 2872 2836                                                                
CONECT 2873 2836                                                                
CONECT 2874 2836                                                                
CONECT 2875 2848                                                                
CONECT 2876 2839                                                                
CONECT 2877 2838    
"""

    rmol = MolFromPDBBlock(block, sanitize=False)
    amol = pdb_dict_to_pdb(pdb_string_to_pdb_dict(block)).model().molecule()

    result = remove_nonpdb_bonds(rmol, amol)

    SanitizeMol(result)
