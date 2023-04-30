#!/usr/bin/env python3
"""Generic infinite network type."""

from __future__ import annotations

from logging import debug
from math import comb
from typing import Any

from gudhi.simplex_tree import SimplexTree
import networkx as nx
import numpy as np
import numpy.typing as npt
from tqdm import tqdm

from distribution.empirical_distribution import EmpiricalDistribution
from network.network import Network
from network.property import BaseNetworkProperty, DerivedNetworkProperty


class FiniteNetwork(Network):
    """Base class representing a data set or a simulated finite network."""

    def __init__(self, max_dimension: int) -> None:
        """Construct an empty network."""
        super().__init__(max_dimension)
        self._is_persistence_computed: bool = False

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        simplicial_complex = SimplexTree()
        for node in tqdm(self.graph.nodes, desc='Inserting nodes from graph to simplical complex', delay=10):
            simplicial_complex.insert((node,))
        for edge in tqdm(self.graph.edges, desc='Inserting edges from graph to simplical complex', delay=10):
            simplicial_complex.insert(edge)

        self._interactions = self.graph.edges
        self._facets = self.graph.edges
        self.simplicial_complex = simplicial_complex

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

    def get_component(self, component_index: int) -> FiniteNetwork:
        """Reduce the network to the specified component only."""
        reduced_network = FiniteNetwork(self.max_dimension)

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
        debug(f'Calculating {property_type.name}...')

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
        elif property_type == BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1:
            property_value = self._calculate_higher_order_degree_distribution(1)
        elif property_type == BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2:
            property_value = self._calculate_higher_order_degree_distribution(2)
        elif property_type == BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_3:
            property_value = self._calculate_higher_order_degree_distribution(3)
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

        debug('finished.')

        return property_value

    def _calculate_average_degree(self) -> float:
        num_of_nodes = self.graph.number_of_nodes()
        num_of_edges = self.graph.number_of_edges()
        average_degree = num_of_edges / num_of_nodes
        return average_degree

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

    def _calculate_higher_order_degree_distribution(
        self,
        simplex_dimension: int,
        neighbor_dimension: int | None = None,
    ) -> EmpiricalDistribution:

        if simplex_dimension >= self.max_dimension:
            # simplex_dimension is at most max_dimension - 1 as higher dimensions make no sense
            return EmpiricalDistribution([])

        neighbor_dimension = simplex_dimension + 1 if neighbor_dimension is None else neighbor_dimension
        assert neighbor_dimension > simplex_dimension

        if simplex_dimension == 0 and neighbor_dimension == 1:
            # ordinry degrees are faster to extract from the networkx package
            simplex_degrees = self.graph.degree()
            degree_sequence = sorted((degree for _, degree in simplex_degrees), reverse=True)
        else:
            degree_sequence = self._calc_degree_sequence(simplex_dimension, neighbor_dimension)

        return EmpiricalDistribution(degree_sequence)

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

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simlex dimension {simplex_dimension}.'

        # extract facets with dimension higher or equal to neighbor dimension
        facets = filter(lambda facet: len(facet) - 1 >= neighbor_dimension, self.facets)

        degree_sequence: list[int] = []
        for facet in tqdm(
            facets,
            desc=f'Calculating degree sequence ({simplex_dimension} - {neighbor_dimension})...',
            delay=10
        ):
            facet_dimension = len(facet) - 1
            num_of_simplices_in_facet = comb(facet_dimension + 1, simplex_dimension + 1)
            num_of_neighbors_in_facet = comb(
                facet_dimension - simplex_dimension,
                neighbor_dimension - simplex_dimension
            )
            degree_sequence += [num_of_neighbors_in_facet] * num_of_simplices_in_facet

        return degree_sequence

    @property
    def simplices(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        simplices_with_filtration = self.simplicial_complex.get_simplices()
        simplices = [simplex for simplex, _ in simplices_with_filtration]
        return simplices

    @property
    def is_persistence_computed(self) -> bool:
        """Return if persistence is computed on the simplicial complex."""
        return self._is_persistence_computed
