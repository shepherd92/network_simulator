"""Helper functions for logging."""

from logging import debug, getLogger


def log_function_name(func: callable):
    """Log decorated function name.

    Args:
        func (_type_): called function
    """
    def wrap(*args, **kwargs):

        debug(f"Calling {func.__name__}")
        result = func(*args, **kwargs)
        debug(f"Finished {func.__name__}")

        logger = getLogger()
        logger.handlers[0].flush()

        return result

    return wrap
