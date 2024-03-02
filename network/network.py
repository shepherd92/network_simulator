#!/usr/bin/env python3
"""Generic network type."""

from __future__ import annotations

from itertools import chain
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

# pylint: disable-next=no-name-in-module
from cpp_modules.build.network import Network as CppNetwork
from distribution.empirical_distribution import EmpiricalDistribution


class Network:
    """Base class representing a data set or a simulated network."""

    def __init__(self, max_dimension: int) -> None:
        """Construct an empty network."""
        assert isinstance(max_dimension, int)
        self._simplex_dimension_distribution: EmpiricalDistribution | None = None
        self._cpp_network = CppNetwork(max_dimension)
        self._graph: nx.Graph() | None = None
        self._digraph: nx.DiGraph() | None = None
        self._vertex_positions: dict[int, tuple[float, ...]] | None = None
        self._interaction_positions: dict[int, tuple[float, ...]] | None = None

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        raise NotImplementedError

    def generate_clique_complex_from_graph(self) -> None:
        """Given graph, generate the clique complex and store it in the simplicial complex member."""
        self.generate_simplicial_complex_from_graph()
        self.expand()

    def expand(self) -> None:
        """Expand simplicial complex to a clique complex."""
        self._cpp_network.expand()

    def filter_to_graph(self) -> None:
        """Filter out those simplices from the simplicial complex that are not present in the graph."""
        graph_nodes_set = list(self.graph.nodes)
        self._cpp_network.keep_only_vertices(graph_nodes_set)
        self.filter_interactions_from_graph()

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
        graph.add_nodes_from([
            simplex[0] for simplex in self._cpp_network.get_skeleton(0)
        ])

        graph.add_weighted_edges_from([
            [*simplex, weight]
            for simplex, weight in self._cpp_network.get_skeleton(1)
            if len(simplex) == 2
        ])
        return graph

    def add_simplices(self, simplices: list[list[int]]) -> None:
        """Insert simplex to the simplicial complex."""
        self._cpp_network.add_simplices(simplices)

    def get_simplices_by_dimension(self, dimension: int) -> list[list[int]]:
        """Return all simplices by their dimension."""
        result: list[list[int]] = [
            list(sorted(simplex))
            for simplex in tqdm(self.simplices, desc='Extracting simplices by dimension', delay=10)
            if len(simplex) - 1 == dimension
        ]
        return result

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
        self._interactions = None
        self.graph = None
        self._digraph = None

    def _calculate_simplex_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the estimated number of simplices for each dimension.

        Facets might have nonempty intersections which are estimated to be empty.
        """
        return EmpiricalDistribution(self._cpp_network.get_simplex_dimension_distribution())

    def _calculate_facet_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of facets for each dimension."""
        facet_dimension_distribution = self._calc_dimension_distribution(self.facets)
        return facet_dimension_distribution

    def _calculate_vertex_interaction_degree_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each vertex."""
        vertices_in_interactions = np.array(list(chain(*self.interactions)))
        degree_sequence = np.unique(vertices_in_interactions, return_counts=True)[1]

        num_of_zero_degree_vertices = self.num_vertices - len(degree_sequence)
        interaction_degree_distribution = EmpiricalDistribution(
            list(degree_sequence) + [0] * num_of_zero_degree_vertices)
        return interaction_degree_distribution

    def _calculate_interaction_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each dimension."""
        interaction_dimension_distribution = self._calc_dimension_distribution(self.interactions)
        return interaction_dimension_distribution

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simplex dimension {simplex_dimension}.'

        assert self._cpp_network.num_vertices() > 0, 'Simplicial complex is empty.'

        degree_sequence: list[int] = self._cpp_network.calc_degree_sequence(simplex_dimension, neighbor_dimension)

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

    @staticmethod
    def _calc_dimension_distribution(simplices: list[list[int]]) -> EmpiricalDistribution:
        dimensions = [len(simplex) - 1 for simplex in simplices]
        distribution = EmpiricalDistribution(dimensions)
        return distribution

    @property
    def simplices(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        raise NotImplementedError

    @property
    def num_simplices(self) -> int:
        """Return the number of simplices in the simplicial complex."""
        return self._cpp_network.num_simplices()

    @property
    def num_vertices(self) -> int:
        """Return the number of nodes in the graph."""
        return self._cpp_network.num_vertices()

    @property
    def max_dimension(self) -> int:
        """Get the maximum dimension of the network."""
        return self._cpp_network.max_dimension()

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
    def interactions(self) -> list[list[int]]:
        """Get network interactions."""
        return self._interactions

    @interactions.setter
    def interactions(self, value: list[list[int]]) -> None:
        """Setter for interactions."""
        self._interactions = value

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
