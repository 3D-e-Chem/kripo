import pkg_resources

from rdkit.Chem import AllChem, MolToSmiles, SanitizeMol


class Reactor:
    """Reactor, facilitator of reactions

    Attributes:
        steps: Maximum number of times a product is reacted

    """
    def __init__(self, steps=99):
        self.reactions = []
        self.steps = steps
        self.load_reactions()

    def load_reactions(self):
        stream = pkg_resources.resource_stream('kripo', 'data/combined.smirks')
        for line in stream:
            cols = line.decode('ascii').strip().split()
            if len(cols) > 0 and cols[0]:
                smirk = cols[0]
                reaction = AllChem.ReactionFromSmarts(smirk)
                self.reactions.append(reaction)

    def react(self, reactant):
        products = set()
        product_smiles = set()
        n = self.steps
        new_mols = [reactant]
        while n > 0 and new_mols != []:
            mols = new_mols
            new_mols = []
            for mol in mols:
                SanitizeMol(mol)
                for reaction in self.reactions:
                    for ps in reaction.RunReactants((mol,)):
                        q = ps[0]
                        SanitizeMol(q)
                        smile = MolToSmiles(q)
                        if smile not in product_smiles:
                            new_mols.append(q)
                            product_smiles.add(smile)
                            products.add(q)

            n -= 1

        return products
