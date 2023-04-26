"""Association Rules: Basic
"""

from romancal.associations import Association
from romancal.associations.lib.constraint import ConstraintTrue
from romancal.associations.registry import RegistryMarker


@RegistryMarker.rule
class Rule_1(Association):
    """Basic rule"""

    def __init__(self, version_id=None):
        self.constraints = ConstraintTrue()
        super().__init__(version_id=version_id)
        self.data["members"] = []

    def _add(self, item):
        self.data["members"].append(item)

    def is_item_member(self, item):
        """Check if item is already a member of this association

        Parameters
        ----------
        item: dict
            The item to add.

        Returns
        -------
        is_item_member: bool
            True if item is a member.
        """
        return item in self.data["members"]
