"""Functionality required by DMS outside of core functionality of the package

Example: Certain tests are required by DMS to log in a specific way. The
decorator for this is defined in this module.
"""

from stpipe import log as stpipe_log


def log_result(requirement, message, result):
    """Log individual test results that relate to a requirement

    Parameters
    ----------
    requirement : str
        The requirement being logged. I.e "DMS363"

    message : str
        Message describing what is being tested

    result : bool
        The result of the test
    """
    logger = stpipe_log.delegator.log
    result_text = "PASS" if result else "FAIL"
    log_msg = f"{requirement} MSG: {message}.......{result_text}"
    logger.info(log_msg)
