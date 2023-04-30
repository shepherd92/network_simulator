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

    def _build_simplicial_complex(self) -> None:
        """Build a simplicial complex based on the loaded data."""
        simplices = [
            [3, 4, 5, 6],
            [2, 3, 4],
            [1, 2], [2, 3], [1, 7], [2, 7],
            [0]
        ]
        # pylint: disable-next=attribute-defined-outside-init
        self._interactions = simplices
        self._facets = []
        self.add_simplices(simplices)

    def _build_graph(self) -> None:
        """Build a simple networkx graph."""
        self.generate_graph_from_simplicial_complex()

    def __str__(self) -> str:
        """Return a string representation based on the data set properties."""
        return '\n'.join([
            'Name: Test',
        ])
