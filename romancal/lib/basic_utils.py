"""General utility objects"""

import numpy as np
from roman_datamodels.datamodels import AssociationsModel
from roman_datamodels.dqflags import group, pixel


def bytes2human(n):
    """Convert bytes to human-readable format

    Taken from the `psutil` library which references
    http://code.activestate.com/recipes/578019

    Parameters
    ----------
    n : int
        Number to convert

    Returns
    -------
    readable : str
        A string with units attached.

    Examples
    --------
    >>> bytes2human(10000)
        '9.8K'

    >>> bytes2human(100001221)
        '95.4M'
    """
    symbols = ("K", "M", "G", "T", "P", "E", "Z", "Y")
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return f"{value:.1f}{s}"
    return f"{n}B"


def is_fully_saturated(model):
    """
    Check to see if all data pixels are flagged as saturated.
    """

    if np.all(np.bitwise_and(model.groupdq, group.SATURATED) == group.SATURATED):
        return True
    elif np.all(np.bitwise_and(model.pixeldq, pixel.SATURATED) == pixel.SATURATED):
        return True

    return False


def is_association(asn_data):
    """
    Test if an object is an association by checking for required fields
    """
    return AssociationsModel.is_association(asn_data)


class LoggingContext:
    """Logging context manager
    Keep logging configuration within a context
    Based on the Python 3 Logging Cookbook example
    Parameters
    ==========
    logger: logging.Logger
        The logger to modify.
    level: int
        The log level to set.
    handler: logging.Handler
        The handler to use.
    close: bool
        Close the handler when done.
    """

    def __init__(self, logger, level=None, handler=None, close=True):
        self.logger = logger
        self.level = level
        self.handler = handler
        self.close = close

        self.old_level = None

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)
        if self.handler:
            self.logger.addHandler(self.handler)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)
        if self.handler:
            self.logger.removeHandler(self.handler)
        if self.handler and self.close:
            self.handler.close()
        # implicit return of None => don't swallow exceptions


def recarray_to_ndarray(x, to_dtype="<f8"):
    """
    Convert a structured array to a 2D ndarray.

    Parameters
    ----------
    x : np.record
        Structured array
    to_dtype : str
        Cast all columns in `x` to this dtype in the output ndarray.

    Returns
    -------
    array : np.ndarray
        Numpy array (without labeled columns).
    """
    names = x.dtype.names
    astype = [(name, to_dtype) for name in names]
    return np.asarray(x.astype(astype).view(to_dtype).reshape((-1, len(names))))


def parse_visitID(visit_id):
    """Utility to parse the visit_id into its components

    Input:
       visit_id as a string

    Output:
      program number
      execution plan number
      pass number
      segment number
      observation number
      visit number
    """

    visit_id_parts = dict(
        [
            ("Program", visit_id[0:5]),
            ("Execution", visit_id[5:7]),
            ("Pass", visit_id[7:10]),
            ("Segment", visit_id[10:13]),
            ("Observation", visit_id[13:16]),
            ("Visit", visit_id[16:20]),
        ]
    )

    return visit_id_parts
