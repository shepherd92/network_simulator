#!/usr/bin/env python3
"""Generic network type."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from logging import debug
from math import comb
from typing import Any

from gudhi.simplex_tree import SimplexTree
import networkx as nx
import numpy as np
import numpy.typing as npt

from distribution.empirical_distribution import EmpiricalDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty


class Network:
    """Base class representing a data set or a simulated network."""

    def __init__(self, max_dimension: int) -> None:
        """Construct an empty network."""
        self._interactions: list[set[int]] = []
        self._simplicial_complex = SimplexTree()
        self._is_persistence_computed: bool = False
        self._graph = nx.Graph()
        self._digraph = nx.DiGraph()
        self._max_dimension = max_dimension

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        simplicial_complex = SimplexTree()
        for node in self.graph.nodes:
            simplicial_complex.insert((node,))
        for edge in self.graph.edges:
            simplicial_complex.insert(edge)

        self.interactions = self.graph.edges
        self.simplicial_complex = simplicial_complex

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
        if len(simplex) == 0:
            return

        self._is_persistence_computed = False

        skeleton = self._get_simplex_skeleton_for_max_dimension(simplex)
        for face in skeleton:
            self.simplicial_complex.insert(face, filtration)

    def get_component(self, component_index: int) -> Network:
        """Reduce the network to the specified component only."""
        reduced_network = Network(self.max_dimension)

        if component_index != -1:
            components = sorted(nx.connected_components(self._graph), key=len, reverse=True)
            reduced_network.graph = self._graph.subgraph(components[component_index])
            reduced_network.generate_clique_complex_from_graph()

        return reduced_network

    def reduce_to_component(self, component_index: int) -> None:
        """Reduce the network to the specified component only."""
        if component_index != -1:
            components = sorted(nx.connected_components(self._graph), key=len, reverse=True)
            self._graph = self._graph.subgraph(components[component_index])
            self.filter_simplicial_complex_from_graph()

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

    def calc_network_summary(
        self,
        properties_to_calculate: list[BaseNetworkProperty.Type]
    ) -> dict[BaseNetworkProperty.Type, Any]:
        """Calculate the summary of the network."""
        summary: dict[BaseNetworkProperty.Type, Any] = {
            property_type: self.calc_base_property(property_type)
            for property_type in properties_to_calculate
        }
        return summary

    def calc_scalar_property(
        self,
        scalar_property_params: DerivedNetworkProperty,
    ) -> float | int:
        """Calculate scalar network properties.

        Note: the parameter of this function contains a callable calculator method.
        This makes this function impossible to use in parallelized settings.
        """
        base_network_property_type = scalar_property_params.source_base_property.property_type
        source_base_property = self.calc_base_property(base_network_property_type)
        scalar_property_value = scalar_property_params.calculator(source_base_property)
        return scalar_property_value

    def calc_base_property(self, property_type: BaseNetworkProperty.Type) -> Any:
        """Return a base property of the network."""
        if property_type == BaseNetworkProperty.Type.NUM_OF_NODES:
            property_value = self.graph.number_of_nodes()
        elif property_type == BaseNetworkProperty.Type.NUM_OF_EDGES:
            property_value = self.graph.number_of_edges()
        elif property_type == BaseNetworkProperty.Type.AVERAGE_DEGREE:
            property_value = self._calculate_average_degree()
        elif property_type == BaseNetworkProperty.Type.MAX_DEGREE:
            oridnary_degree_distributions = self._calculate_degree_distribution()
            property_value = oridnary_degree_distributions.domain.max_
        elif property_type == BaseNetworkProperty.Type.DEGREE_DISTRIBUTION:
            property_value = self._calculate_degree_distribution()
        elif property_type == BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION:
            property_value = self._calculate_in_degree_distribution()
        elif property_type == BaseNetworkProperty.Type.OUT_DEGREE_DISTRIBUTION:
            property_value = self._calculate_out_degree_distribution()
        elif property_type == BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTIONS:
            higher_order_degree_distributions = self._calculate_all_higher_order_degree_distributions()
            property_value = higher_order_degree_distributions
        elif property_type == BaseNetworkProperty.Type.AVG_CLUSTERING:
            property_value = nx.average_clustering(self.graph)
        elif property_type == BaseNetworkProperty.Type.NUM_OF_CONNECTED_COMPONENTS:
            property_value = nx.number_connected_components(self.graph)
        elif property_type == BaseNetworkProperty.Type.DIMENSION:
            property_value = self.simplicial_complex.dimension()
        elif property_type == BaseNetworkProperty.Type.NUM_OF_SIMPLICES:
            property_value = self.simplicial_complex.num_simplices()
        elif property_type == BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_interaction_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.SIMPLEX_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_simplex_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.FACET_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_facet_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.BETTI_NUMBERS:
            property_value = self._calculate_betti_numbers()
        elif property_type == BaseNetworkProperty.Type.PERSISTENCE:
            property_value = tuple(self.simplicial_complex.persistence(persistence_dim_max=True))
        else:
            raise NotImplementedError(
                f'Requested property type {property_type} is not available.'
            )

        return property_value

    def _calculate_average_degree(self) -> float:
        num_of_nodes = self.graph.number_of_nodes()
        num_of_edges = self.graph.number_of_edges()
        average_degree = num_of_edges / num_of_nodes
        return average_degree

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

    def _compute_persistence(self) -> None:
        """Compute persistence of the simplicial complex."""
        debug('Computing persistence...')
        # self.simplicial_complex.compute_persistence(persistence_dim_max=self.max_dimension)
        self.simplicial_complex.compute_persistence()
        self._is_persistence_computed = True
        debug('Persistence computed.')

    def _calculate_degree_distribution(self) -> EmpiricalDistribution:
        degree_sequence = [degree for _, degree in self.graph.degree()]
        return EmpiricalDistribution(degree_sequence)

    def _calculate_in_degree_distribution(self) -> EmpiricalDistribution:
        degree_sequence = [degree for _, degree in self.digraph.in_degree()]
        return EmpiricalDistribution(degree_sequence)

    def _calculate_out_degree_distribution(self) -> EmpiricalDistribution:
        degree_sequence = [degree for _, degree in self.digraph.out_degree()]
        return EmpiricalDistribution(degree_sequence)

    def _calculate_all_higher_order_degree_distributions(self) -> dict[int, EmpiricalDistribution]:

        result: dict[int, EmpiricalDistribution] = {}
        for simplex_dimension in range(1, self.max_dimension):
            # the for loop goes up to max_dimension - 1 as higher dimensions make no sense
            ho_distribution = self._calculate_higher_order_degree_distribution(simplex_dimension)
            result[simplex_dimension] = ho_distribution

        return result

    def _calculate_higher_order_degree_distribution(
        self,
        simplex_dimension: int,
        neighbor_dimension: int | None = None,
    ) -> EmpiricalDistribution:

        neighbor_dimension = simplex_dimension + 1 if neighbor_dimension is None else neighbor_dimension
        assert neighbor_dimension > simplex_dimension

        if simplex_dimension == 0 and neighbor_dimension == 1:
            # ordinry degrees are faster to extract from the networkx package
            simplex_degrees = self.graph.degree()
            degree_sequence = sorted((degree for _, degree in simplex_degrees), reverse=True)
        else:
            degree_sequence = self._calc_degree_sequence(simplex_dimension, neighbor_dimension)

        return EmpiricalDistribution(degree_sequence)

    def _calc_degree_sequence(
        self,
        simplex_dimension: int,
        neighbor_dimension: int
    ) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simlex dimension {simplex_dimension}.'

        # extract facets with dimension higher or equal to neighbor dimension
        facets = filter(lambda facet: len(facet) > neighbor_dimension, self.extract_facets())

        degree_sequence: list[int] = []
        for facet in facets:
            facet_dimension = len(facet) - 1
            num_of_simplices_in_facet = comb(facet_dimension + 1, simplex_dimension + 1)
            num_of_neighbors_in_facet = comb(
                facet_dimension - simplex_dimension,
                neighbor_dimension - simplex_dimension
            )
            degree_sequence += [num_of_neighbors_in_facet] * num_of_simplices_in_facet

        return degree_sequence

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

    def _calculate_betti_numbers(self) -> npt.NDArray[np.int_]:
        """Calculate the Betti numbers for different dimensions."""
        if not self.is_persistence_computed:
            self._compute_persistence()

        if self.simplicial_complex.dimension() != 0:
            betti_numbers = self.simplicial_complex.betti_numbers()
        else:
            betti_numbers = [self.graph.number_of_nodes()]

        betti_number_array = np.zeros((self.max_dimension, 2))
        betti_number_array[:, 0] = range(self.max_dimension)

        for dimension, betti_number in enumerate(betti_numbers):
            betti_number_array[dimension, 1] = betti_number

        return betti_number_array

    @staticmethod
    def _calc_dimension_distribution(simplices: list[list[int]]) -> EmpiricalDistribution:
        dimensions = [len(simplex) - 1 for simplex in simplices]
        distribution = EmpiricalDistribution(dimensions)
        return distribution

    @property
    def interactions(self) -> list[set[int]]:
        """Get the interactions associated to the data set."""
        return self._interactions

    @interactions.setter
    def interactions(self, value: list[set[int]]) -> None:
        self._interactions = value

    @interactions.deleter
    def interactions(self) -> None:
        del self._interactions

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
    def simplices(self) -> list[set[int]]:
        """Get the simplices associated to the network."""
        simplices_with_filtration = self.simplicial_complex.get_simplices()
        simplices = [simplex for simplex, _ in simplices_with_filtration]
        return simplices

    @property
    def dimension(self) -> int:
        """Get the dimension of the simplicial complex."""
        return self.simplicial_complex.dimension

    @property
    def max_dimension(self) -> int:
        """Get the maximum dimension of the network."""
        return self._max_dimension

    @property
    def is_persistence_computed(self) -> bool:
        """Return if persistence is computed on the simplicial complex."""
        return self._is_persistence_computed

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
