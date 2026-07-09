"""Load an Association from a file or object"""

from ._association import _Association


def load_asn(
    serialized,
    validate=True,
):
    """Load an Association from a file or object

    Parameters
    ----------
    serialized : object
        The serialized form of the association.

    validate : bool
        Validate against the class' defined schema, if any.

    Returns
    -------
    dict
    The association data

    Raises
    ------
    AssociationNotValidError
        Cannot create or validate the association.

    Notes
    -----
    The `serialized` object can be either a string or
    a file object containing the string.
    """
    return _Association.load(serialized, validate=validate)
