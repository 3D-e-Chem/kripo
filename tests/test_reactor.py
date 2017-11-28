import pytest

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
            'c1(ccc(cc1)*)*',
            'c1c(ccc(c1C(F)(F)F)Cl)*',
            'c1c(ccc(c1C(F)(F)F)Cl)NC(=O)N*',
            'c1c(ccc(c1C(F)(F)F)Cl)NC(=O)Nc1ccc(cc1)*',
            'c1c(ccc(c1C(F)(F)F)Cl)NC(=O)Nc1ccc(cc1)O*',
            'c1c(ccc(c1C(F)(F)F)Cl)NC(=O)Nc1ccc(cc1)Oc1ccnc(c1)*',
            'c1c(ccc(c1C(F)(F)F)Cl)NC(=O)Nc1ccc(cc1)Oc1ccnc(c1)C(=O)NC',
            'c1c(ccc(c1*)Cl)*',
            'c1c(ccc(c1*)Cl)NC(=O)N*',
            'c1c(ccc(c1*)Cl)NC(=O)Nc1ccc(cc1)*',
            'c1c(ccc(c1*)Cl)NC(=O)Nc1ccc(cc1)O*',
            'c1c(ccc(c1*)Cl)NC(=O)Nc1ccc(cc1)Oc1ccnc(c1)*',
            'c1c(ccc(c1*)Cl)NC(=O)Nc1ccc(cc1)Oc1ccnc(c1)C(=O)NC',
            'c1(ccc(cc1)O*)*',
            'c1(ccc(cc1)Oc1ccnc(c1)*)*',
            'c1(ccc(cc1)Oc1ccnc(c1)C(=O)NC)*',
            '*c1ccnc(c1)*',
            '*c1ccnc(c1)C(=O)NC',
            '*C(F)(F)F',
            '*C(=O)NC',
            '*NC(=O)N*',
            '*NC(=O)Nc1ccc(cc1)*',
            '*NC(=O)Nc1ccc(cc1)O*',
            '*NC(=O)Nc1ccc(cc1)Oc1ccnc(c1)*',
            '*NC(=O)Nc1ccc(cc1)Oc1ccnc(c1)C(=O)NC',
            '*O*',
            '*Oc1ccnc(c1)*',
            '*Oc1ccnc(c1)C(=O)NC',
        },
        id='BAX'
    ),
])
def test_react(reactant, expected_products):
    reactor = Reactor()
    reactant_mol = MolFromSmiles(reactant)
    products = reactor.react(reactant_mol)

    products = mols2smiles(products)
    assert products == expected_products
