#!/usr/bin/env python
"""This module specifies helper functions."""


def static_vars(**kwargs):
    """Specify static variables (decorator)."""
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate
