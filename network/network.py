#!/usr/bin/env python3
"""Generic network type."""

from __future__ import annotations

from itertools import chain
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd

# pylint: disable=no-name-in-module
from cpp_modules.build.network import FiniteNetwork as CppFiniteNetwork
from cpp_modules.build.network import InfiniteNetwork as CppInfiniteNetwork
# pylint: enable=no-name-in-module
from distribution.empirical_distribution import EmpiricalDistribution


class Network:
    """Base class representing a data set or a simulated network."""

    def __init__(self) -> None:
        """Construct an empty network."""
        self._cpp_network = None
        self._graph: nx.Graph | None = None
        self._digraph: nx.DiGraph | None = None
        self._vertex_positions: dict[int, tuple[float, ...]] | None = None
        self._interaction_positions: dict[int, tuple[float, ...]] | None = None

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        cpp_network_type = type(self._cpp_network)
        self._cpp_network = cpp_network_type(self.max_dimension, self.graph.nodes, self.graph.edges)

    def generate_clique_complex_from_graph(self) -> None:
        """Given graph, generate the clique complex and store it in the simplicial complex member."""
        self.generate_simplicial_complex_from_graph()
        self.expand()

    def create_simplicial_complex(self) -> None:
        """Create a simplicial complex."""
        self._cpp_network.create_simplicial_complex()

    def expand(self) -> None:
        """Expand simplicial complex to a clique complex."""
        self._cpp_network.expand(self.max_dimension)

    def filter_to_graph(self) -> None:
        """Filter out those simplices from the simplicial complex that are not present in the graph."""
        self._cpp_network.keep_only_vertices(list(self.graph.nodes))

    def filter_interactions_from_graph(self) -> None:
        """Filter out interactions that are not present in the graph."""
        graph_nodes_set = set(self.graph.nodes)
        filtered_interactions: list[list[int]] = [
            list(set(interaction) & graph_nodes_set)
            for interaction in self.interactions
        ]
        # filter out empty lists
        self.interactions = list(filter(None, filtered_interactions))

    def _generate_graph_from_simplicial_complex(self) -> nx.Graph:
        """Set graph to represent the simplicial complex."""
        # assert self.simplicial_complex.num_vertices() != 0, 'Simplicial complex is not built.'

        graph = nx.Graph()
        graph.add_nodes_from(self.vertices)
        graph.add_edges_from([
            [*simplex]
            for simplex in self._cpp_network.get_skeleton(1)
            if len(simplex) == 2
        ])
        return graph

    def num_of_vertices_in_component(self, component_index: int) -> int:
        """Return the number of vertices in a component."""
        components = sorted(nx.connected_components(self.graph), key=len, reverse=True)
        if component_index >= len(components):
            return 0

        return self.graph.subgraph(components[component_index]).number_of_nodes()

    def get_info_as_dict(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        raise NotImplementedError

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        info = self.get_info_as_dict()
        data_frame = pd.DataFrame(info, index=[0])
        data_frame.to_csv(save_path, index=False)

    def _reset(self) -> None:
        self._cpp_network.reset()
        self.graph = None
        self._digraph = None

    def _calculate_simplex_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the estimated number of simplices for each dimension.

        Facets might have nonempty intersections which are estimated to be empty.
        """
        return EmpiricalDistribution(self.cpp_network.calc_simplex_dimension_distribution())

    def _calculate_facet_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of facets for each dimension."""
        return EmpiricalDistribution(self.cpp_network.calc_facet_dimension_distribution())

    def _calculate_vertex_interaction_degree_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each vertex."""
        vertices_in_interactions = np.array(list(chain(*self.interactions)))
        degree_sequence = np.unique(vertices_in_interactions, return_counts=True)[1]

        num_of_zero_degree_vertices = self.num_simplices(0) - len(degree_sequence)
        interaction_degree_distribution = EmpiricalDistribution(
            list(degree_sequence) + [0] * num_of_zero_degree_vertices)
        return interaction_degree_distribution

    def _calculate_interaction_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each dimension."""
        return EmpiricalDistribution(self.cpp_network.calc_interaction_dimension_distribution())

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simplex dimension {simplex_dimension}.'

        degree_sequence: list[int] = self.cpp_network.calc_degree_sequence(simplex_dimension, neighbor_dimension)

        # degree_sequence: list[int] = [
        #     len(self.simplicial_complex.get_cofaces(simplex, neighbor_dimension - simplex_dimension))
        #     for simplex in tqdm(self.simplices, desc='Calculating degree sequence', delay=10)
        #     if len(simplex) == simplex_dimension + 1
        # ]
        # degree_sequence: list[int] = calc_degree_sequence(
        #     self.facets,
        #     simplex_dimension,
        #     neighbor_dimension
        # )
        return degree_sequence

    @property
    def simplices(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        return self._cpp_network.get_simplices()

    @property
    def cpp_network(self) -> CppFiniteNetwork | CppInfiniteNetwork | None:
        """Return the cpp network."""
        return self._cpp_network

    def num_simplices(self, dimension: int = -1) -> int:
        """Return the number of simplices in the simplicial complex."""
        if dimension == -1:
            return self._cpp_network.num_simplices()
        elif dimension == 0:
            return self._cpp_network.num_vertices()
        elif dimension == 1:
            return self._num_edges()
        else:
            simplex_dimension_value_counts = self._calculate_simplex_dimension_distribution().calc_value_counts()
            return 0 \
                if len(simplex_dimension_value_counts[simplex_dimension_value_counts[:, 0] == dimension]) == 0 \
                else simplex_dimension_value_counts[simplex_dimension_value_counts[:, 0] == dimension][0, 1]

    def _num_edges(self) -> int:
        """Return the number of edges in the graph."""
        raise NotImplementedError

    @property
    def max_dimension(self) -> int:
        """Get the maximum dimension of the network."""
        return self._cpp_network.max_dimension

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
    def facets(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        return self._cpp_network.get_facets()

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
