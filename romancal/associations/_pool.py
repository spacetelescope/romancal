"""
Association Pools
"""

from collections import UserDict


class PoolRow(UserDict):
    """A row from an AssociationPool

    Class to create a copy of an AssociationPool row without copying
    all of the astropy.Table.Row private attributes.
    """

    def __init__(self, init=None):
        dict_init = dict(init)
        super().__init__(dict_init)
        try:
            self.meta = init.meta
        except AttributeError:
            self.meta = dict()
