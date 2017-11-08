import pkg_resources

from rdkit.Chem import AllChem
from rdkit.Chem import Recap, Mol


class Reactor:
    """Reactor, facilitator of reactions"""
    # def __init__(self):
    #     self.reactions = []
    #     self.load_reactions()
    #
    # def load_reactions(self):
    #     stream = pkg_resources.resource_stream('kripo', 'data/combined.smirks')
    #     for line in stream:
    #         cols = line.decode('ascii').strip().split()
    #         if len(cols) > 0 and cols[0]:
    #             smirk = cols[0]
    #             # TODO fix smirks so they are OK for rdkit
    #             reaction = AllChem.ReactionFromSmarts(smirk)
    #             self.reactions.append(reaction)
    #
    # def react(self, reactant):
    #     products = []
    #     for reaction in self.reactions:
    #         products.extend(reaction.RunReactants((reactant, )))
    #     # TODO make products unique
    #     return products

    def react(self, reactant: Mol):
        """Reacts reactant to products

        Args:
            reactant (Mol): Reactant

        Returns:
            List[Mol]: Products

        """
        hierarch = Recap.RecapDecompose(reactant)
        return [c.mol for c in hierarch.GetAllChildren().values()]

