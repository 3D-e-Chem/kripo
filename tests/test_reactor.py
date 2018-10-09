import pytest
from kripo.pdb import pdb_from_file, ligands
from rdkit.Chem import AllChem, MolFromMolBlock, MolToMolBlock

from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles

from kripo.reactor import Reactor, embed_r_groups


def mols2smiles(mols):
    return {MolToSmiles(m) for m in mols}


@pytest.mark.parametrize("reactant, expected_products", [
    pytest.param(
        'C',
        set(),
        id='C'
    ),
    pytest.param(
        'CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(c3)C(F)(F)F)cc2)ccn1',
        {
            '*NC(=O)Nc1ccc(Oc2ccnc(C(=O)NC)c2)cc1',
            '*Oc1ccc(NC(=O)Nc2ccc(Cl)c(*)c2)cc1',
            '*c1ccc(Cl)c(*)c1',
            '*Oc1ccc(NC(=O)Nc2ccc(Cl)c(C(F)(F)F)c2)cc1',
            '*NC(=O)Nc1ccc(O*)cc1',
            '*C(=O)NC',
            '*NC(=O)N*',
            '*NC(=O)Nc1ccc(Cl)c(*)c1',
            '*Oc1ccc(*)cc1',
            '*c1ccc(NC(=O)Nc2ccc(Cl)c(C(F)(F)F)c2)cc1',
            '*NC(=O)Nc1ccc(*)cc1',
            '*c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(*)c3)cc2)ccn1',
            '*NC(=O)Nc1ccc(Cl)c(C(F)(F)F)c1',
            '*c1ccc(Oc2ccnc(*)c2)cc1',
            '*c1ccnc(C(=O)NC)c1',
            '*c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1',
            '*c1ccc(*)cc1',
            '*c1cc(NC(=O)Nc2ccc(Oc3ccnc(C(=O)NC)c3)cc2)ccc1Cl',
            '*c1ccc(Oc2ccnc(C(=O)NC)c2)cc1',
            '*C(F)(F)F',
            '*c1ccnc(*)c1',
            '*c1ccc(Cl)c(C(F)(F)F)c1',
            '*O*',
            '*Oc1ccnc(C(=O)NC)c1',
            '*Oc1ccnc(*)c1',
            '*NC(=O)Nc1ccc(Oc2ccnc(*)c2)cc1',
            '*c1ccc(NC(=O)Nc2ccc(Cl)c(*)c2)cc1'
        },
        id='BAX'
    ),
])
def test_react(reactant, expected_products):
    reactor = Reactor()
    reactant_mol = MolFromSmiles(reactant)
    AllChem.EmbedMolecule(reactant_mol, AllChem.ETKDG())
    products = reactor.react(reactant_mol)

    products = mols2smiles(products)
    assert products == expected_products


@pytest.fixture
def bax_mol():
    return MolFromMolBlock('''
RDKit          3D

 48 50  0  0  0  0  0  0  0  0999 V2000
   -1.9780    2.2590   25.9480 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2280    3.4280   25.8570 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6980    4.5910   26.4520 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9200    4.6600   27.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4480    2.7040   24.2850 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3380    2.3780   22.9000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6560    2.2710   21.6650 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3330    1.8220   20.5270 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7070    1.4910   20.6340 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3820    1.5940   21.8510 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6920    2.0440   22.9670 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4640    1.0880   19.5620 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7130   -1.0290   18.8880 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4170   -1.9280   17.8630 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0410   -0.3090   16.2360 C   0  0  1  0  0  0  0  0  0  0  0  0
    4.3280    0.6230   17.2290 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1010   -1.0340   13.8920 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.2530   -1.0040   12.7050 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7100    3.5350   27.1680 C   0  0  1  0  0  0  0  0  0  0  0  0
   -3.2190    2.2650   26.5840 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0500    0.9970   26.6980 C   0  0  1  0  0  0  0  0  0  0  0  0
   -5.2380    1.1420   26.2110 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4800    0.0010   26.1220 F   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2250    0.6490   27.9280 F   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2760    3.5990   28.0420 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.0220    3.4870   25.2640 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.7480    2.8090   24.0510 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2620    1.9460   23.6160 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1550    0.2360   18.5580 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5810   -1.5370   16.5780 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.1310    0.0000   14.7620 C   0  0  1  0  0  0  0  0  0  0  0  0
    4.1550    1.1730   14.3910 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1120    5.4400   26.4100 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3160    3.2190   24.7430 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5870    4.1890   25.5880 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0780   -2.8780   18.0820 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.6590    1.5710   16.9900 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.6630   -1.8240   14.0630 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6580    2.5250   21.6040 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2350    5.5450   27.5270 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.7140   -0.0800   12.6750 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5610   -1.8190   12.7400 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8620   -1.0910   11.8290 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.1900    2.1340   23.8660 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.3800    1.3390   21.9200 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8410    1.7330   19.6240 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6140    1.3830   25.5420 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.6030   -1.3050   19.8760 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1 20  1  0
  1 47  1  0
  2  3  1  0
  2 26  1  0
  3  4  1  0
  3 33  1  0
  4 19  1  0
  4 40  1  0
  5 26  1  0
  5 27  1  0
  5 28  1  0
  6  7  1  0
  6 11  1  0
  6 27  1  0
  7  8  1  0
  7 39  1  0
  8  9  1  0
  8 46  1  0
  9 10  1  0
  9 12  1  0
 10 11  1  0
 10 45  1  0
 11 44  1  0
 12 29  1  0
 13 14  1  0
 13 29  1  0
 13 48  1  0
 14 30  1  0
 14 36  1  0
 15 16  1  0
 15 30  1  1
 15 31  1  0
 16 29  1  0
 16 37  1  0
 17 18  1  0
 17 31  1  0
 17 38  1  0
 18 41  1  0
 18 42  1  0
 18 43  1  0
 19 20  1  0
 19 25  1  1
 20 21  1  0
 21 22  1  6
 21 23  1  0
 21 24  1  0
 26 35  1  0
 27 34  1  0
 31 32  1  6
M  END
    ''')


def test_embed_r_groups__ROR(bax_mol):
    fragment = MolFromMolBlock('''
RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 R   0  0  0  0  0  1  0  0  0  0  0  0
    4.4640    1.0880   19.5620 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 R   0  0  0  0  0  1  0  0  0  0  0  0
  1  2  1  0
  3  2  1  0
M  END
    ''')
    embed_r_groups(fragment, bax_mol)

    expected = '''
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    3.7070    1.4910   20.6340 R   0  0  0  0  0  1  0  0  0  0  0  0
    4.4640    1.0880   19.5620 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1550    0.2360   18.5580 R   0  0  0  0  0  1  0  0  0  0  0  0
  1  2  1  0
  3  2  1  0
M  END
'''
    assert MolToMolBlock(fragment) == expected


def test_embed_r_groups__noparentatom(bax_mol):
    fragment = MolFromMolBlock('''
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 R   0  0  0  0  0  1  0  0  0  0  0  0
    1.1111    1.1111   11.1111 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4480    2.7040   24.2850 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3160    3.2190   24.7430 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0220    3.4870   25.2640 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2620    1.9460   23.6160 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 R   0  0  0  0  0  1  0  0  0  0  0  0
    0.5870    4.1890   25.5880 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  3  2  1  0
  2  4  1  0
  3  5  1  0
  3  6  1  0
  7  5  1  0
  5  8  1  0
M  END''')

    with pytest.raises(LookupError) as e:
        embed_r_groups(fragment, bax_mol)
        assert 'group not found in parent' in str(e.value)

