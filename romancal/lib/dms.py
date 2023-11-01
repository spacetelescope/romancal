"""Functionality required by DMS outside of core functionality of the package

Example: Certain tests are required by DMS to log in a specific way. The
decorator for this is defined in this module.
"""
from functools import wraps

def result_logger(message, logger):
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
                result_text = 'PASS' if result else 'FAIL'
                logger.info(message + result_text)

            try:
                test_function(*args, **kwargs)
            except Exception:
                log_result(False)
                raise
            else:
                log_result(True)

        return inner

    return decorator

