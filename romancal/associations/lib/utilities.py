"""General Utilities"""

import logging
from ast import literal_eval
from functools import wraps

from numpy.ma import masked

from romancal.lib.basic_utils import parse_visitID as parse_visitID

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def return_on_exception(exceptions=(Exception,), default=None):
    """Decorator to force functions raising exceptions to return a value

    Parameters
    ----------
    exceptions: (Exception(,...))
        Tuple of exceptions to catch

    default: obj
        The value to return when a specified exception occurs
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except exceptions as err:
                logger.debug(
                    "Caught exception %s in function %s, forcing return value of %s",
                    err,
                    func,
                    default,
                )
                return default

        return wrapper

    return decorator


def evaluate(value):
    """Evaluate a value

    Parameters
    ----------
    value : str
        The string to evaluate.

    Returns
    -------
    type or str
        The evaluation. If the value cannot be
        evaluated, the value is simply returned
    """
    try:
        evaled = literal_eval(value)
    except (ValueError, SyntaxError):
        evaled = value
    return evaled


def getattr_from_list(adict, attributes, invalid_values=None):
    """Retrieve value from dict using a list of attributes

    Parameters
    ----------
    adict : dict
        dict to retrieve from

    attributes : list
        List of attributes

    invalid_values : set
        A set of values that essentially mean the
        attribute does not exist.

    Returns
    -------
    (attribute, value)
        Returns the value and the attribute from
        which the value was taken.

    Raises
    ------
    KeyError
        None of the attributes are found in the dict.
    """
    if invalid_values is None:
        invalid_values = set()

    for attribute in attributes:
        try:
            result = adict[attribute]
        except KeyError:
            continue
        else:
            if result is masked:
                continue
            if result not in invalid_values:
                return attribute, result
            else:
                continue
    else:
        raise KeyError(f"Object has no attributes in {attributes}")


@return_on_exception(exceptions=(KeyError,), default=None)
def getattr_from_list_nofail(*args, **kwargs):
    """Call getattr_from_list without allows exceptions.

    If the specified exceptions are caught, return `default`
    instead.

    Parameters
    ----------
    See `getattr_from_list`
    """
    return getattr_from_list(*args, **kwargs)


def is_iterable(obj):
    """General iterator check"""
    return (
        not isinstance(obj, str)
        and not isinstance(obj, tuple)
        and hasattr(obj, "__iter__")
    )


def mk_level3_asn_name(
    visit_id, output_file_root, filter_id, release_product, product_type, patch_name
):
    """Construct an association file name based on the visit_id and product"""

    parsed_visit_id = parse_visitID(visit_id)

    sep = "_"

    product_name_mapping = {
        "visit": "v"
        + parsed_visit_id["Execution"]
        + parsed_visit_id["Pass"]
        + parsed_visit_id["Segment"]
        + parsed_visit_id["Observation"]
        + parsed_visit_id["Visit"],
        "daily": "d"
        + parsed_visit_id["Execution"]
        + parsed_visit_id["Pass"]
        + parsed_visit_id["Segment"],
        "pass": "p" + parsed_visit_id["Execution"] + parsed_visit_id["Pass"],
        "full": "full",
        "user": "user",
    }

    pr_name = product_name_mapping.get(product_type, "unknown")

    asn_file_name = (
        output_file_root
        + sep
        + release_product
        + sep
        + pr_name
        + sep
        + patch_name
        + sep
        + filter_id
    )

    return asn_file_name
