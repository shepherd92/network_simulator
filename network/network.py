#!/usr/bin/env python3
"""Generic network type."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any

from gudhi.simplex_tree import SimplexTree
import networkx as nx
import numpy as np
import numpy.typing as npt
import pandas as pd
from tqdm import tqdm

from cpp_modules.build.simplicial_complex import extract_facets
from distribution.empirical_distribution import EmpiricalDistribution


class Network:
    """Base class representing a data set or a simulated network."""

    def __init__(self, max_dimension: int) -> None:
        """Construct an empty network."""
        assert isinstance(max_dimension, int)
        self._max_dimension = max_dimension

        self._interactions: list[list[int]] = []
        self._facets: list[list[int]] = []
        self._simplex_dimension_distribution: EmpiricalDistribution | None = None
        self._simplicial_complex = SimplexTree()
        self._graph: nx.Graph() | None = None
        self._digraph: nx.DiGraph() | None = None

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        raise NotImplementedError

    def generate_clique_complex_from_graph(self) -> None:
        """Given graph, generate the clique complex and store it in the simplicial complex member."""
        self.generate_simplicial_complex_from_graph()
        self.expand()

    def expand(self) -> None:
        """Expand simplicial complex to a clique complex."""
        self.simplicial_complex.expansion(self._max_dimension)

    def filter_simplicial_complex_from_graph(self) -> None:
        """Filter out those simplices from the simplicial complex that are not present in the graph."""
        simplicial_complex = SimplexTree()
        for simplex in tqdm(self.simplices, desc='Filtering simplices in the graph', delay=10):
            if all(node in self.graph.nodes for node in simplex):
                # if the first node is in the graph, the whole simplex is there as well
                simplicial_complex.insert(simplex)
        self.simplicial_complex = simplicial_complex

        graph_nodes_set = set(self.graph.nodes)
        interactions: list[list[int]] = [
            list[set(interaction) & graph_nodes_set]
            for interaction in self.interactions
        ]
        # filter out empty lists
        self.interactions = list(filter(None, interactions))

    def _generate_graph_from_simplicial_complex(self) -> nx.Graph:
        """Set graph to represent the simplicial complex."""
        # assert self.simplicial_complex.num_vertices() != 0, 'Simplicial complex is not built.'

        graph = nx.Graph()
        graph.add_nodes_from([
            simplex[0] for simplex, weight in self.simplicial_complex.get_skeleton(0)
        ])

        graph.add_weighted_edges_from([
            [*simplex, weight]
            for simplex, weight in self.simplicial_complex.get_skeleton(1)
            if len(simplex) == 2
        ])
        return graph

    def add_simplices_batch(self, simplices: npt.NDArray[np.int]) -> None:
        """Add simplices to the simplicial complex in batch."""
        assert simplices.shape[1] - 1 <= self.max_dimension, f'Simplices have too high dimension {simplices.shape[1]}.'
        self.simplicial_complex.insert_batch(simplices.T, np.zeros((simplices.shape[0],)))

    def add_simplices(self, simplices: list[list[int]]) -> None:
        """Insert a simplex to the simplicial complex."""
        for simplex in tqdm(simplices, desc='Adding simplices', delay=10):
            self.add_simplex(simplex)

    def add_simplex(self, simplex: list[int]) -> None:
        """Insert a simplex to the simplicial complex.

        Add the skeleton of the simplex as its dimension is too high.
        """
        raise NotImplementedError

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
        components = sorted(nx.connected_components(self._graph), key=len, reverse=True)
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
        self._interactions = []
        self._facets = []
        self._simplex_dimension_distribution = None
        self._simplicial_complex = SimplexTree()
        self._graph = None
        self._digraph = None

    def _extract_facets(self) -> list[list[int]]:
        """Return the facets of the simplicial complex."""
        if self._interactions:
            print(f'Extracting facets from {len(self._interactions)} interactions.')
            facets: list[list[int]] = extract_facets(self._interactions)
        else:
            facets = list(nx.find_cliques(self.graph))

        return facets

    def _calculate_simplex_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the estimated number of simplices for each dimension as a list.

        Facets might have nonempty intersections which are estimated to be empty.
        """
        simplex_dimension_distribution = self._calc_dimension_distribution(self.simplices)
        return simplex_dimension_distribution

    def _calculate_facet_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of facets for each dimension as a list."""
        facet_dimension_distribution = self._calc_dimension_distribution(self.facets)
        return facet_dimension_distribution

    def _calculate_interaction_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each dimension as a list."""
        interaction_dimension_distribution = self._calc_dimension_distribution(self.interactions)
        return interaction_dimension_distribution

    def _get_simplex_skeleton_for_max_dimension(self, simplex: list[int]) -> npt.NDArray[np.int_]:
        simplex_dimension = len(simplex) - 1
        if simplex_dimension <= self.max_dimension:
            skeleton = np.array([simplex])
        else:
            skeleton = np.array(list(combinations(simplex, self.max_dimension + 1)))

        return skeleton

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:
        raise NotImplementedError

    def _calc_degrees(
        self,
        simplex_dimension: int,
        neighbor_dimension: int
    ) -> list[tuple[set[int], int]]:

        higher_order_degrees: list[tuple[set[int], int]] = []
        for simplex in tqdm(self.simplices, desc='Calculating higher-order degrees', delay=10):
            if len(simplex) - 1 == simplex_dimension:
                cofaces = [
                    face
                    for face, _ in self.simplicial_complex.get_cofaces(simplex, simplex_dimension)
                ]
                degree = list(map(len, cofaces)).count(neighbor_dimension + 1)
                higher_order_degrees.append((simplex, degree))

        return higher_order_degrees

    @staticmethod
    def _calc_dimension_distribution(simplices: list[list[int]]) -> EmpiricalDistribution:
        dimensions = [len(simplex) - 1 for simplex in simplices]
        distribution = EmpiricalDistribution(dimensions)
        return distribution

    @property
    def simplicial_complex(self) -> SimplexTree:
        """Get the simplicial complex associated to the network."""
        return self._simplicial_complex

    @simplicial_complex.setter
    def simplicial_complex(self, value: SimplexTree) -> None:
        self._simplicial_complex = value
        self._interactions = None
        self._facets = None
        self._simplex_dimension_distribution = None
        self._is_persistence_computed = False
        self._num_of_edges = None

    @simplicial_complex.deleter
    def simplicial_complex(self) -> None:
        del self._simplicial_complex

    @property
    def simplices(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        raise NotImplementedError

    @property
    def num_simplices(self) -> int:
        """Return the number of simplices in the simplicial complex."""
        return self.simplicial_complex.num_simplices()

    @property
    def num_vertices(self) -> int:
        """Return the number of nodes in the graph."""
        return self.simplicial_complex.num_vertices()

    @property
    def dimension(self) -> int:
        """Get the dimension of the simplicial complex."""
        return self.simplicial_complex.dimension

    @property
    def max_dimension(self) -> int:
        """Get the maximum dimension of the network."""
        return self._max_dimension

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
        if not self._facets:
            self._facets = self._extract_facets()
        return self._facets

    @property
    def interactions(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        if not self._interactions:
            return self.facets
        return self._interactions

    @property
    def simplex_dimension_distribution(self) -> EmpiricalDistribution:
        """Get the simplices associated to the network."""
        if not self._simplex_dimension_distribution:
            self._simplex_dimension_distribution = self._calculate_simplex_dimension_distribution()
        return self._simplex_dimension_distribution
