#!/usr/bin/env python3
"""Base class for all data sets."""

from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import NewType

# pylint: disable-next=no-name-in-module
from cpp_plugin.build.release.cpp_plugin import FiniteHypergraph as CppFiniteHypergraph
from network.finite_hypergraph import FiniteHypergraph
from tools.logging_helper import log_function_name


DerivedNetworkProperty = NewType('DerivedNetworkProperty', None)


class Dataset(FiniteHypergraph):
    """Base class representing a data set."""

    class Type(Enum):
        """Represent the available data set types."""

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
        weighted: bool

    @log_function_name
    def __init__(self, data_set_properties: Parameters) -> None:
        """Create data set without loading data."""
        self._data_set_properties = data_set_properties
        self._read_data()
        vertices = self._get_vertices_from_data()
        interactions = self._get_interactions_from_data()
        cpp_network = CppFiniteHypergraph(
            data_set_properties.max_dimension,
            vertices, interactions,
            data_set_properties.weighted
        )

        super().__init__(cpp_network)
        self.reduce_to_component(self._data_set_properties.component_index_from_largest)

    @log_function_name
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

        source_base_property = self.calc_base_property(scalar_property_params.source_base_property)
        scalar_property_value = calculator(source_base_property)
        return scalar_property_value

    @log_function_name
    def _read_data(self) -> None:
        """Read data from the disk."""
        raise NotImplementedError

    @log_function_name
    def _get_vertices_from_data(self) -> list[int]:
        """Return vertices of the loaded data."""
        raise NotImplementedError

    @log_function_name
    def _get_interactions_from_data(self) -> list[list[int]]:
        """Get interactions from the data."""
        raise NotImplementedError
