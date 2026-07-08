"""Reprocessing List"""

from enum import Enum

__all__ = [
    "ListCategory",
    "ProcessList",
]


class ListCategory(Enum):
    """Categories in the list"""

    RULES = 0
    BOTH = 1
    EXISTING = 2
    NONSCIENCE = 3


class ProcessList:
    """A Process list

    Parameters
    ----------
    items : [item[, ...]]
        The list of items to process

    rules : [Association[, ...]]
        List of rules to process the items against.

    work_over : int
        What the reprocessing should work on:
        - `ProcessList.EXISTING`:   Only existing associations
        - `ProcessList.RULES`:      Only on the rules to create new associations
        - `ProcessList.BOTH`:       Compare to both existing and rules
        - `ProcessList.NONSCIENCE`: Only on non-science items

    only_on_match : bool
        Only use this object if the overall condition
        is True.

    trigger_constraints : [Constraint[,...]]
        The constraints that created the ProcessList

    trigger_rules : [Association[,...]]
        The association rules that created the ProcessList
    """

    _str_attrs = (
        "rules",
        "work_over",
        "only_on_match",
        "trigger_constraints",
        "trigger_rules",
    )

    def __init__(
        self,
        items=None,
        rules=None,
        work_over=ListCategory.BOTH,
        only_on_match=False,
        trigger_constraints=None,
        trigger_rules=None,
    ):
        self.items = items
        self.rules = rules
        self.work_over = work_over
        self.only_on_match = only_on_match
        self.trigger_constraints = (
            set(trigger_constraints) if trigger_constraints else set()
        )
        self.trigger_rules = set(trigger_rules) if trigger_rules else set()

    @property
    def hash(self):
        """Create a unique hash"""
        return (tuple(self.rules), self.work_over, self.only_on_match)

    def update(self, process_list, full=False):
        """Update with information from ProcessList

        Attributes from `process_list` are added to self's attributes. If `not
        full`, the attributes `rules`, 'work_over`, and `only_on_match` are not
        taken.

        Note that if `full`, destructive action will occur with respect to
        `work_over` and `only_on_match`.

        Parameters
        ----------
        process_list : ProcessList
            The source process list to absorb.

        full : bool
            Include the hash attributes `rules`, `work_over`, and `only_on_match`.
        """
        self.items += process_list.items
        self.trigger_constraints.update(process_list.trigger_constraints)
        self.trigger_rules.update(process_list.trigger_rules)
        if full:
            self.rules += process_list.rules
            self.work_over = process_list.work_over
            self.only_on_match = process_list.only_on_match

    def __str__(self):
        result = f"{self.__class__.__name__}(n_items: {len(self.items)}, { ({str_attr: getattr(self, str_attr) for str_attr in self._str_attrs}) })"
        return result
