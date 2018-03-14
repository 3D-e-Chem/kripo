from typing import List, Tuple

from ..fragment import Fragment
from ..pdb import MAX_CONTACT_DISTANCE
from .vector import distance_between_positions


class Feature:
    """Pharmacophore feature

    Attributes:
        kind (str): Kind of feature. Must one of kripodb.pharmacophores.FEATURE_TYPE_KEYS
        position (Tuple(float, float, float): X, Y, Z position of feature

    """
    def __init__(self, kind, position):
        self.kind = kind
        self.position = position

    def distance(self, other):
        """Euclidean distance of feature position to other.

        Args:
            other (Feature): Other feature

        Returns:
            float: Distance in angstrom.
        """
        pos1 = self.position
        pos2 = other.position
        return distance_between_positions(pos1, pos2)

    def in_contact_with(self, fragment_positions: List[Tuple[float, float, float]]):
        """Test whether feature is closer than MAX_CONTACT_DISTANCE to fragment

        Args:
            fragment_positions : The fragment atom positions

        Returns:
            bool: True if feature is in contact with fragment.
        """
        for a in fragment_positions:
            dist = distance_between_positions(self.position, a)
            if dist < MAX_CONTACT_DISTANCE:
                return True
        return False

    def __repr__(self) -> str:
        p = self.position
        return 'Feature(\'{0}\', [{1:.4}, {2:.4}, {3:.4}])'.format(self.kind, p[0], p[1], p[2])

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other) -> bool:
        return repr(self) == repr(other)
