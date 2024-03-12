#!/usr/bin/env python3
"""Test data set."""

from dataclasses import dataclass

from data_set.data_set import DataSet


class TestDataSet(DataSet):
    """This class represents the "Bianconi Data Set"."""

    @dataclass
    class Parameters(DataSet.Parameters):
        """Properties to load the test data set."""

    def _read_data(self) -> None:
        """Load data from the disk for further processing."""
        return

    def _get_interactions(self) -> list[list[int]]:
        """Interactions of the loaded data."""
        return [
            # [0, 1, 2, 3], [1, 2, 3, 4],
            [3, 4, 5, 6],
            [2, 3, 4], [4, 5, 8], [4, 6, 8], [5, 6, 8],
            [1, 2], [1, 7], [2, 7],
            [0]
        ]

    def __str__(self) -> str:
        """Return a string representation based on the data set properties."""
        return '\n'.join([
            'Name: Test',
        ])
