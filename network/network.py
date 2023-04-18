#!/usr/bin/env python3
"""Generic network type."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations

from gudhi.simplex_tree import SimplexTree
import networkx as nx

from distribution.empirical_distribution import EmpiricalDistribution


class Network:
    """Base class representing a data set or a simulated network."""

    def __init__(self, max_dimension: int) -> None:
        """Construct an empty network."""
        self._interactions: list[set[int]] = []
        self._simplicial_complex = SimplexTree()
        self._graph = nx.Graph()
        self._digraph = nx.DiGraph()
        assert isinstance(max_dimension, int)
        self._max_dimension = max_dimension

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        raise NotImplementedError

    def generate_clique_complex_from_graph(self) -> None:
        """Given graph, generate the clique complex and store it in the simplicial complex member."""
        cliques = list(nx.find_cliques(self.graph))
        self.simplicial_complex = SimplexTree()
        self.add_simplices(cliques)
        self.interactions = [set(clique) for clique in cliques]

    def filter_simplicial_complex_from_graph(self) -> None:
        """Filter those simplices from the simplicial complex that are not present in the graph."""
        simplicial_complex = SimplexTree()
        for simplex in self.simplices:
            if all(node in self.graph.nodes for node in simplex):
                # if the first node is in the graph, the whole simplex is there as well
                simplicial_complex.insert(simplex)
        self.simplicial_complex = simplicial_complex

        graph_nodes_set = set(self.graph.nodes)
        interactions: list[set[int]] = [
            interaction & graph_nodes_set
            for interaction in self.interactions
        ]
        # filter out empty lists
        self.interactions = list(filter(None, interactions))

    def generate_graph_from_simplicial_complex(self) -> None:
        """Set graph to represent the simplicial complex."""
        assert self.simplicial_complex.num_vertices() != 0, 'Simplicial complex is not built.'

        for simplex, weight in self.simplicial_complex.get_skeleton(1):
            if len(simplex) == 1:
                self.graph.add_node(*simplex)

        self.graph.add_weighted_edges_from([
            [*simplex, weight]
            for simplex, weight in self.simplicial_complex.get_skeleton(1)
            if len(simplex) == 2
        ])

        self.digraph = self.graph.to_directed()

    def add_simplices(self, simplices: list[list[int]]) -> None:
        """Insert a simplex to the simplicial complex."""
        for simplex in simplices:
            self.add_simplex(simplex)

    def add_simplex(self, simplex: list[int], filtration: float = 0.) -> None:
        """Insert a simplex to the simplicial complex.

        Add the skeleton of the simplex as its dimension is too high.
        """
        raise NotImplementedError

    def extract_facets(self) -> list[set[int]]:
        """Return the facets of the simplicial complex."""
        facets = list(filter(
            lambda interaction: not any(
                interaction < other_interaction
                for other_interaction in self.interactions
            ), self.interactions
        ))

        return facets

    def get_simplices_by_dimension(self) -> dict[int, list[list[int]]]:
        """Return all simplices by their dimension."""
        result: dict[int, list[list[int]]] = defaultdict(list)
        for simplex in self.simplices:
            dimension = len(simplex) - 1
            result[dimension].append(list(sorted(simplex)))

        return result

    def _calculate_simplex_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the estimated number of simplices for each dimension as a list.

        Facets might have nonempty intersections which are estimated to be empty.
        """
        simplex_dimension_distribution = self._calc_dimension_distribution(self.simplices)
        return simplex_dimension_distribution

    def _calculate_facet_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of facets for each dimension as a list."""
        facets = self.extract_facets()
        facet_dimension_distribution = self._calc_dimension_distribution(facets)
        return facet_dimension_distribution

    def _calculate_interaction_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each dimension as a list."""
        interaction_dimension_distribution = self._calc_dimension_distribution(self.interactions)
        return interaction_dimension_distribution

    def _get_simplex_skeleton_for_max_dimension(self, simplex: list[int]) -> list[list[int]]:
        simplex_dimension = len(simplex) - 1
        if simplex_dimension <= self.max_dimension:
            skeleton = [simplex]
        else:
            skeleton = list(combinations(simplex, self.max_dimension + 1))

        return skeleton

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:
        raise NotImplementedError

    def _calc_degrees(
        self,
        simplex_dimension: int,
        neighbor_dimension: int
    ) -> list[tuple[set[int], int]]:

        higher_order_degrees: list[tuple[set[int], int]] = []
        for simplex in self.simplices:
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
        self._is_persistence_computed = False

    @simplicial_complex.deleter
    def simplicial_complex(self) -> None:
        del self._simplicial_complex

    @property
    def simplices(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        raise NotImplementedError

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
        return self._graph

    @graph.setter
    def graph(self, value: nx.Graph) -> None:
        self._graph = value

    @property
    def digraph(self) -> nx.DiGraph:
        """Get the simple directed graph associated to the network."""
        return self._digraph

    @digraph.setter
    def digraph(self, value: nx.DiGraph) -> None:
        self._digraph = value
