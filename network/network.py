#!/usr/bin/env python3
"""Generic network type."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd

# pylint: disable-next=no-name-in-module
from cpp_plugin.build.release.cpp_plugin import (
    FiniteCliqueComplex as CppFiniteCliqueComplex,
    InfiniteCliqueComplex as CppInfiniteCliqueComplex
)
from distribution.empirical_distribution import EmpiricalDistribution
from network.property import BaseNetworkProperty
from tools.logging_helper import log_function_name


class Network:
    """Base class representing a data set or a simulated network."""

    def __init__(self) -> None:
        """Construct an empty network."""
        self._cpp_network = None
        self._graph: nx.Graph | None = None
        self._digraph: nx.DiGraph | None = None
        self._vertex_positions: dict[int, tuple[float, ...]] | None = None
        self._interaction_positions: dict[int, tuple[float, ...]] | None = None

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

    @log_function_name
    def num_of_vertices_in_component(self, component_index: int) -> int:
        """Return the number of vertices in a component."""
        components = sorted(nx.connected_components(self.graph), key=len, reverse=True)
        if component_index >= len(components):
            return 0

        return self.graph.subgraph(components[component_index]).number_of_nodes()

    def info(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        raise NotImplementedError

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        info = self.info()
        data_frame = pd.DataFrame(info, index=[0])
        data_frame.to_csv(save_path, index=False)

    @log_function_name
    def _calculate_simplex_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the estimated number of simplices for each dimension."""
        return EmpiricalDistribution(self.cpp_network.calc_simplex_dimension_distribution())

    def _calculate_vertex_interaction_degree_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each vertex."""
        raise NotImplementedError

    def _calculate_interaction_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each dimension."""
        raise NotImplementedError

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simplex dimension {simplex_dimension}.'

        return self.cpp_network.calc_coface_degree_sequence(simplex_dimension, neighbor_dimension)

    @property
    def cpp_network(self) -> CppFiniteCliqueComplex | CppInfiniteCliqueComplex | None:
        """Return the cpp network."""
        return self._cpp_network

    @cpp_network.setter
    def cpp_network(self, value: CppFiniteCliqueComplex | CppInfiniteCliqueComplex) -> None:
        """Set the cpp network."""
        self._cpp_network = value

    @log_function_name
    def num_simplices(self, dimension: int) -> int:
        """Return the number of simplices in the simplicial complex."""
        return self.cpp_network.num_simplices(dimension)

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
    def digraph(self) -> nx.DiGraph:
        """Get the simple directed graph associated to the network."""
        if self._digraph is None:
            self._digraph = self.graph.to_directed()
        return self._digraph

    @property
    def vertices(self) -> list[int]:
        """Get network vertices."""
        return self._cpp_network.vertices

    @vertices.setter
    def vertices(self, value: list[int]) -> None:
        """Setter for vertices."""
        self._cpp_network.vertices = value

    @property
    def interactions(self) -> list[list[int]]:
        """Get network interactions."""
        return self._cpp_network.interactions

    @interactions.setter
    def interactions(self, value: list[list[int]]) -> None:
        """Setter for interactions."""
        self._cpp_network.interactions = value

    @property
    def vertex_positions(self) -> dict[int, tuple[float, ...]]:
        """Return vertex positions to plot."""
        return self._vertex_positions

    @vertex_positions.setter
    def vertex_positions(self, value: dict[int, tuple[float, ...]]) -> None:
        """Setter of vertex positions."""
        self._vertex_positions = value

    @property
    def interaction_positions(self) -> dict[int, tuple[float, ...]]:
        """Return interaction positions to plot."""
        return self._interaction_positions

    @interaction_positions.setter
    def interaction_positions(self, value: dict[int, tuple[float, ...]]) -> None:
        """Setter of interaction positions."""
        self._interaction_positions = value

    def __str__(self) -> str:
        """Return a string representation based on the network properties."""
        return '\n'.join([
            f'{key}: {item}'
            for key, item in self.info().items()
        ])
