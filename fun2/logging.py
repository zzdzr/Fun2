"""
This module provides common decorators for logging and skipping execution
based on the existence of a file or directory.
"""

import os
import functools
import logging

def log_step(func):
    """Decorator that logs the start and end of a function execution.

    Args:
        func (Callable): The function to wrap.

    Returns:
        Callable: The wrapped function.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logging.info("Starting %s...", func.__name__)
        try:
            result = func(*args, **kwargs)
            logging.info("Finished %s successfully.", func.__name__)
            return result
        except Exception as e:
            logging.error("Error in %s: %s", func.__name__, e, exc_info=True)
            raise
    return wrapper

def skip_if_exists(get_path):
    """Decorator that skips function execution if the specified path exists.

    Args:
        get_path (Callable): A function that accepts the instance (usually `self`)
            and returns the file or directory path to check.

    Returns:
        Callable: The wrapped function that is skipped if the path exists.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            path = get_path(args[0])
            if os.path.exists(path):
                logging.info("Path '%s' already exists. Skipping %s.", path, func.__name__)
                return path
            return func(*args, **kwargs)
        return wrapper
    return decorator