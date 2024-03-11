#!/usr/bin/env python3
"""Base class for all data sets."""

from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path

from network.finite_network import FiniteNetwork
from network.property import DerivedNetworkProperty


class DataSet(FiniteNetwork):
    """Base class representing a data set."""

    class Type(Enum):
        """Represent the available data set types."""

        BIANCONI = auto()
        NATURE = auto()
        ARXIV = auto()
        TEST = auto()
        INVALID = auto()

    @dataclass
    class Parameters:
        """Represent the necessary properties to load a dataset."""

        location: Path
        max_dimension: int
        max_simplex_dimension: int
        component_index_from_largest: int

    def __init__(self, data_set_properties: Parameters) -> None:
        """Create data set without loading data."""
        self._data_set_properties = data_set_properties
        self._read_data()
        vertices = self._get_vertices()
        interactions = self._get_interactions()

        super().__init__(data_set_properties.max_dimension, vertices, interactions)
        self.reduce_to_component(self._data_set_properties.component_index_from_largest)

    def calc_scalar_property(
        self,
        scalar_property_params: DerivedNetworkProperty,
    ) -> float | int:
        """Calculate scalar network properties.

        Note: the parameter of this function contains a callable calculator method.
        This makes this function impossible to use in parallelized settings.
        """
        calculator = scalar_property_params.calculator_data_set \
            if scalar_property_params.calculator_data_set is not None \
            else scalar_property_params.calculator_default

        base_network_property_type = scalar_property_params.source_base_property.property_type
        source_base_property = self.calc_base_property(base_network_property_type)
        scalar_property_value = calculator(source_base_property)
        return scalar_property_value

    def _read_data(self) -> None:
        """Read data from the disk."""
        raise NotImplementedError

    def _get_vertices(self) -> list[int]:
        """Build a simplicial complex based on the loaded data."""
        raise NotImplementedError

    def _get_interactions(self) -> list[list[int]]:
        """Get interactions."""
        raise NotImplementedError
