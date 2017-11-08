from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from atomium.files.pdbdict2pdb import pdb_dict_to_pdb
from atomium.structures.chains import Site
from rdkit.Chem import MolFromPDBBlock

from kripo.site import protonate, to_pdb_block, protonate_rdkit_molecule, protonate_site_from_rdkit, chain_of_site


def test_protonate(site_3heg_bax: Site):
    site = protonate(site_3heg_bax)

    hydrogens = site.atoms(element='H')
    assert len(hydrogens) > 0


def test_to_pdb_block(site_3heg_bax: Site):
    block = to_pdb_block(site_3heg_bax)

    assert 'ATOM' in block


def test_protonate_site_from_rdkit():
    block = """ATOM     15  N   ALA A  51      -1.046   1.075  -0.025  1.00 39.82           N  
ATOM     16  CA  ALA A  51      -0.147   0.049  -0.470  1.00 38.51           C  
ATOM     17  C   ALA A  51      -0.640  -1.305  -0.093  1.00 37.91           C  
ATOM     18  O   ALA A  51      -1.694  -1.386   0.519  1.00 38.13           O  
ATOM     19  CB  ALA A  51       1.212   0.225   0.149  1.00 38.89           C  
"""
    pdb = pdb_dict_to_pdb(pdb_string_to_pdb_dict(block))
    # Treat chain as site as they both inherit from AtomicStructure and ResidueStructure
    unprotenated_site = pdb.model().chain()

    unprotonated_mol = MolFromPDBBlock(block)
    protonated_mol = protonate_rdkit_molecule(unprotonated_mol)

    protonated_site = protonate_site_from_rdkit(unprotenated_site, protonated_mol)

    hydrogens = protonated_site.atoms(element='H')
    assert len(hydrogens) == 7


def test_chain_of_site(site_3heg_bax: Site):
    chain = chain_of_site(site_3heg_bax)

    assert chain == 'A'
