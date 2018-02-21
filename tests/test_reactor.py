import pytest
from rdkit.Chem import AllChem

from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles

from kripo.reactor import Reactor


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
            '[*]NC(=O)Nc1ccc(Oc2ccnc(C(=O)NC)c2)cc1',
            '[*]Oc1ccc(NC(=O)Nc2ccc(Cl)c([*])c2)cc1',
            '[*]c1ccc(Cl)c([*])c1',
            '[*]Oc1ccc(NC(=O)Nc2ccc(Cl)c(C(F)(F)F)c2)cc1',
            '[*]NC(=O)Nc1ccc(O[*])cc1',
            '[*]C(=O)NC',
            '[*]NC(=O)N[*]',
            '[*]NC(=O)Nc1ccc(Cl)c([*])c1',
            '[*]Oc1ccc([*])cc1',
            '[*]c1ccc(NC(=O)Nc2ccc(Cl)c(C(F)(F)F)c2)cc1',
            '[*]NC(=O)Nc1ccc([*])cc1',
            '[*]c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c([*])c3)cc2)ccn1',
            '[*]NC(=O)Nc1ccc(Cl)c(C(F)(F)F)c1',
            '[*]c1ccc(Oc2ccnc([*])c2)cc1',
            '[*]c1ccnc(C(=O)NC)c1',
            '[*]c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1',
            '[*]c1ccc([*])cc1',
            '[*]c1cc(NC(=O)Nc2ccc(Oc3ccnc(C(=O)NC)c3)cc2)ccc1Cl',
            '[*]c1ccc(Oc2ccnc(C(=O)NC)c2)cc1',
            '[*]C(F)(F)F',
            '[*]c1ccnc([*])c1',
            '[*]c1ccc(Cl)c(C(F)(F)F)c1',
            '[*]O[*]',
            '[*]Oc1ccnc(C(=O)NC)c1',
            '[*]Oc1ccnc([*])c1',
            '[*]NC(=O)Nc1ccc(Oc2ccnc([*])c2)cc1',
            '[*]c1ccc(NC(=O)Nc2ccc(Cl)c([*])c2)cc1'
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

