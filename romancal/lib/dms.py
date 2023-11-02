"""Functionality required by DMS outside of core functionality of the package

Example: Certain tests are required by DMS to log in a specific way. The
decorator for this is defined in this module.
"""
from functools import wraps
import logging
from metrics_logger.decorators import metrics_logger

from stpipe import log as stpipe_log

def result_logger(requirement, message):
    """Decorator to log assertion status

    Parameters
    ----------
    message : str
        The message to use. Note that the phrase of either "PASS" or "FAIL"
        will be appended to the end of the message.
    """
    def decorator(test_function):

        @wraps(test_function)
        def inner(*args, **kwargs):

            def log_result(result):
                logger = stpipe_log.delegator.log
                result_text = 'PASS' if result else 'FAIL'
                log_msg = f'{requirement} MSG: {message}.......{result_text}'
                logger.info(log_msg)

            metrics_func = metrics_logger(requirement)(test_function)
            try:
                metrics_func(*args, **kwargs)
            except Exception:
                log_result(False)
                raise
            else:
                log_result(True)

        return inner

    return decorator

