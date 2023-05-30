#!/usr/bin/env python3
"""This module is responsible for implement a conditional decorator."""


def conditional_decorator(decorator, condition):
    """Decorate a function if a condition is fulfilled."""
    def decorated_function(inner_function):
        if not condition:
            # Return the function unchanged, not decorated.
            return inner_function
        return decorator(inner_function)
    return decorated_function
