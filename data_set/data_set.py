#!/usr/bin/env python3
"""Base class for all data sets."""

from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path

from network.network import Network


class DataSet(Network):
    """Base class representing a data set."""

    class Type(Enum):
        """Represent the available data set types."""

        BIANCONI: int = auto()
        NATURE: int = auto()
        ARXIV: int = auto()
        TEST: int = auto()
        INVALID: int = auto()

    @dataclass
    class Parameters:
        """Represent the necessary properties to load a dataset."""

        location: Path
        max_dimension: int
        max_simplex_dimension: int
        component_index_from_largest: int

    def __init__(self, data_set_properties: Parameters) -> None:
        """Create data set without loading data."""
        super().__init__(data_set_properties.max_dimension)
        self._data_set_properties = data_set_properties

    def load_data(self) -> None:
        """Load data from the disk for further processing."""
        self._read_data()
        self._build_simplicial_complex()
        self._build_graph()
        self.reduce_to_component(self._data_set_properties.component_index_from_largest)

    def _read_data(self) -> None:
        """Read data from the disk."""
        raise NotImplementedError

    def _build_simplicial_complex(self) -> None:
        """Build a simplicial complex based on the loaded data."""
        raise NotImplementedError

    def _build_graph(self) -> None:
        """Build a simple networkx graph."""
        raise NotImplementedError

    def __str__(self) -> str:
        """Return a string representation based on the data set properties."""
        raise NotImplementedError
