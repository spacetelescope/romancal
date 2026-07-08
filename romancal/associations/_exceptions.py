__all__ = [
    "AssociationError",
    "AssociationNotValidError",
]


class AssociationError(Exception):
    """Basic errors related to Associations"""


class AssociationNotValidError(AssociationError):
    """Given data structure is not a valid association"""
