#!/usr/bin/env python3
"""This module is responsible for patching a the multiprocessing pool istarmap version."""

from multiprocessing import pool


def istarmap(self: pool.Pool, func, iterable, chunksize=1):
    """Starmap-version of imap."""
    self._check_running()
    if chunksize < 1:
        raise ValueError(f'Chunksize must be 1+, not {chunksize:n}')

    task_batches = pool.Pool._get_tasks(func, iterable, chunksize)
    result = pool.IMapIterator(self)
    self._taskqueue.put((
        self._guarded_task_generation(result._job, pool.starmapstar, task_batches),
        result._set_length
    ))
    return (item for chunk in result for item in chunk)


pool.Pool.istarmap = istarmap
