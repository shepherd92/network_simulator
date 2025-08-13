#!/usr/bin/env python3
"""Generic network type."""

from pathlib import Path
from typing import Any, NewType

import networkx as nx
import pandas as pd

from distribution.empirical_distribution import EmpiricalDistribution
from network.property import BaseNetworkProperty
from tools.logging_helper import log_function_name


CppFiniteCliqueComplex = NewType('CppFiniteCliqueComplex', None)
CppInfiniteCliqueComplex = NewType('CppInfiniteCliqueComplex', None)
CppFiniteHypergraph = NewType('CppFiniteHypergraph', None)
CppInfiniteHypergraph = NewType('CppInfiniteHypergraph', None)


class Network:
    """Base class representing a data set or a simulated network."""

    def __init__(self, cpp_network: CppFiniteCliqueComplex | CppInfiniteCliqueComplex | CppFiniteHypergraph | CppInfiniteHypergraph) -> None:
        """Construct an empty network."""
        self._cpp_network = cpp_network
        self._graph: nx.Graph | None = None
        self._vertex_positions: dict[int, tuple[float, ...]] | None = None

    def calc_base_property(self, property_type: BaseNetworkProperty) -> Any:
        """Generate typical properties of the given type."""
        raise NotImplementedError

    @log_function_name
    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        cpp_network_type = type(self.cpp_network)
        self.cpp_network = cpp_network_type(self.max_dimension, self.graph.nodes, self.graph.edges)

    @log_function_name
    def _generate_graph_from_simplicial_complex(self) -> nx.Graph:
        """Set graph to represent the simplicial complex."""
        raise NotImplementedError

    def info(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        return {
            'num_of_vertices': self.num_simplices(0),
            'max_dimension': self.max_dimension,
        }

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        info = self.info()
        data_frame = pd.DataFrame(info, index=[0])
        data_frame.to_csv(save_path, index=False)

    @log_function_name
    def _calculate_simplex_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the estimated number of simplices for each dimension."""
        return EmpiricalDistribution(self.cpp_network.calc_simplex_dimension_distribution())

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simplex dimension {simplex_dimension}.'

        return self.cpp_network.calc_coface_degree_sequence(simplex_dimension, neighbor_dimension)

    @log_function_name
    def num_simplices(self, dimension: int) -> int:
        """Return the number of simplices in the simplicial complex."""
        return self.cpp_network.num_simplices(dimension)

    @property
    def cpp_network(self) -> CppFiniteCliqueComplex | CppInfiniteCliqueComplex | None:
        """Return the cpp network."""
        return self._cpp_network

    @cpp_network.setter
    def cpp_network(self, value: CppFiniteCliqueComplex | CppInfiniteCliqueComplex) -> None:
        """Set the cpp network."""
        self._cpp_network = value

    @property
    def max_dimension(self) -> int:
        """Get the maximum dimension of the network."""
        return self.cpp_network.max_dimension

    @property
    def graph(self) -> nx.Graph:
        """Get the simple graph associated to the network."""
        if self._graph is None:
            self._graph = self._generate_graph_from_simplicial_complex()
        return self._graph

    @graph.setter
    def graph(self, value: nx.Graph) -> None:
        self._graph = value

    @property
    def vertices(self) -> list[int]:
        """Get network vertices."""
        return self._cpp_network.vertices

    @vertices.setter
    def vertices(self, value: list[int]) -> None:
        """Setter for vertices."""
        self._cpp_network.vertices = value

    @property
    def vertex_positions(self) -> dict[int, tuple[float, ...]]:
        """Return vertex positions to plot."""
        return self._vertex_positions

    @vertex_positions.setter
    def vertex_positions(self, value: dict[int, tuple[float, ...]]) -> None:
        """Setter of vertex positions."""
        self._vertex_positions = value

    def __str__(self) -> str:
        """Return a string representation based on the network properties."""
        return '\n'.join([
            f'{key}: {item}'
            for key, item in self.info().items()
        ])
