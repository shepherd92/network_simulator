#!/usr/bin/env python3
"""Generic infinite network type."""

from __future__ import annotations

from logging import debug, info
from typing import Any

import networkx as nx
import numpy as np
import numpy.typing as npt
from tqdm import tqdm

from distribution.empirical_distribution import EmpiricalDistribution
# pylint: disable=no-name-in-module
from cpp_modules.build.network import filter_simplices
from cpp_modules.build.network import FiniteNetwork as CppFiniteNetwork
# pylint: enable=no-name-in-module
from network.network import Network
from network.property import BaseNetworkProperty, DerivedNetworkProperty
from tools.logging_helper import log_function_name


class FiniteNetwork(Network):
    """Base class representing a data set or a simulated finite network."""

    def __init__(self, max_dimension: int, vertices: list[int], interactions: list[list[int]]) -> None:
        """Construct an empty network."""
        super().__init__()
        assert isinstance(max_dimension, int)
        self._cpp_network = CppFiniteNetwork(max_dimension, vertices, interactions)
        self._components: list[FiniteNetwork] | None = None

    @log_function_name
    def get_component(self, component_index: int) -> FiniteNetwork:
        """Return the network of the specified component."""
        if component_index == -1:
            return self
        vertices_in_component = self._get_vertices_in_components()[component_index]
        interactions_in_component = filter_simplices(self.interactions, vertices_in_component)
        return FiniteNetwork(self.max_dimension, vertices_in_component, interactions_in_component)

    @log_function_name
    def reduce_to_component(self, component_index: int) -> None:
        """Reduce the network to the specified component only."""
        if component_index != -1:
            component = self.get_component(component_index)
            if self.vertex_positions is not None:
                self.vertex_positions = {
                    id: position
                    for id, position in self.vertex_positions.items()
                    if id in component.vertices
                }

    @log_function_name
    def calc_network_summary(
        self,
        properties_to_calculate: list[BaseNetworkProperty.Type]
    ) -> dict[BaseNetworkProperty.Type, Any]:
        """Calculate the summary of the network."""
        info('Finite network summary calculation started.')
        summary: dict[BaseNetworkProperty.Type, Any] = {
            property_type: self.calc_base_property(property_type)
            for property_type in properties_to_calculate
        }
        info('Finite network summary calculation finished.')
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
        scalar_property_value = scalar_property_params.calculator_default(source_base_property)
        return scalar_property_value

    @log_function_name
    def collapse(self) -> FiniteNetwork:
        """Collapse edges of the simplicial complex preserving its 1 homology."""
        debug(f'Collapsing network, containing currently: {self.num_simplices(0)} vertices.')
        collapsed_network = FiniteNetwork(self.max_dimension, self.vertices, self.interactions)
        collapsed_network.cpp_network.collapse()
        collapsed_network.expand()
        debug(f'Collapsed containing {collapsed_network.num_simplices(0)} vertices.')
        return collapsed_network

    @log_function_name
    def calc_base_property(self, property_type: BaseNetworkProperty.Type) -> Any:
        """Return a base property of the network."""
        debug(f'Calculating {property_type.name}...')

        if property_type == BaseNetworkProperty.Type.NUM_OF_VERTICES:
            property_value = self.num_simplices(0)
        elif property_type == BaseNetworkProperty.Type.NUM_OF_EDGES:
            property_value = self.num_simplices(1)
        elif property_type == BaseNetworkProperty.Type.NUM_OF_TRIANGLES:
            property_value = self.num_simplices(2)
        elif property_type == BaseNetworkProperty.Type.NUM_OF_INTERACTIONS:
            property_value = len(self.interactions)
        elif property_type == BaseNetworkProperty.Type.EDGES:
            property_value = np.array(self.graph.edges, dtype=int)
        elif property_type == BaseNetworkProperty.Type.AVERAGE_DEGREE:
            property_value = self._calculate_average_degree()
        elif property_type == BaseNetworkProperty.Type.MAX_DEGREE:
            oridnary_degree_distributions = self._calculate_degree_distribution()
            property_value = oridnary_degree_distributions.domain.max_
        elif property_type == BaseNetworkProperty.Type.DEGREE_DISTRIBUTION:
            property_value = self._calculate_degree_distribution()
        elif property_type == BaseNetworkProperty.Type.VERTEX_INTERACTION_DEGREE_DISTRIBUTION:
            property_value = self._calculate_vertex_interaction_degree_distribution()
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
        elif property_type == BaseNetworkProperty.Type.NUM_OF_SIMPLICES:
            property_value = self.num_simplices()
        elif property_type == BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_interaction_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.SIMPLEX_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_simplex_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.FACET_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_facet_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.BETTI_NUMBERS:
            property_value = self._calculate_betti_numbers()
        elif property_type == BaseNetworkProperty.Type.BETTI_NUMBERS_BY_COMPONENT:
            property_value = self._calculate_betti_numbers_in_components()
        elif property_type == BaseNetworkProperty.Type.VERTICES_BY_COMPONENT:
            property_value = self._calculate_vertices_in_components()
        elif property_type == BaseNetworkProperty.Type.PERSISTENCE_PAIRS:
            property_value = self._calc_persistence_pairs()
        else:
            raise NotImplementedError(
                f'Requested property type {property_type} is not available.'
            )

        debug('finished.')

        return property_value

    def get_info_as_dict(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        return {
            'num_of_vertices': self.num_simplices(0),
            'num_of_interactions': len(self.interactions),
            'num_of_simplices': self.num_simplices(),
            'max_dimension': self.max_dimension,
            'num_of_components': len(self._get_vertices_in_components()),
            'num_of_vertices_in_component_0': self.num_of_vertices_in_component(0),
        }

    @log_function_name
    def _calculate_average_degree(self) -> float:
        num_of_nodes = self.graph.number_of_nodes()
        num_of_edges = self.graph.number_of_edges()
        average_degree = num_of_edges / num_of_nodes * 2
        return average_degree

    @log_function_name
    def _calculate_degree_distribution(self) -> EmpiricalDistribution:
        degree_sequence = [degree for _, degree in self.graph.degree()]
        return EmpiricalDistribution(degree_sequence)

    @log_function_name
    def _calculate_in_degree_distribution(self) -> EmpiricalDistribution:
        degree_sequence = [degree for _, degree in self.digraph.in_degree()]
        return EmpiricalDistribution(degree_sequence)

    @log_function_name
    def _calculate_out_degree_distribution(self) -> EmpiricalDistribution:
        degree_sequence = [degree for _, degree in self.digraph.out_degree()]
        return EmpiricalDistribution(degree_sequence)

    @log_function_name
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

    @log_function_name
    def _calculate_betti_numbers_in_components(self) -> npt.NDArray[np.int_]:

        betti_numbers_by_component: list[list[int]] = [
            component.calc_base_property(BaseNetworkProperty.Type.BETTI_NUMBERS)[:, 1]
            for component in tqdm(
                self.components,
                desc='Calculating Betti numbers in components',
                delay=10,
            )
        ]
        return np.array(betti_numbers_by_component)

    @log_function_name
    def _calculate_vertices_in_components(self) -> npt.NDArray[np.int_]:

        vertices_in_components_list: list[int] = [
            len(vertices_in_component)
            for vertices_in_component in tqdm(
                self._get_vertices_in_components(),
                desc='Calculating vertices in components',
                delay=10,
            )
        ]

        vertices_by_component = np.array(vertices_in_components_list)
        return vertices_by_component

    def _calculate_betti_numbers(self) -> list[int]:
        """Calculate the Betti numbers for different dimensions."""
        betti_numbers = self.cpp_network.calc_betti_numbers()
        result = np.c_[
            np.array(range(len(betti_numbers)), dtype=int),
            betti_numbers
        ]
        return result

    @log_function_name
    def _calc_persistence_pairs(self) -> list[tuple[list[int], list[int]]]:
        return self.cpp_network.calc_persistence_pairs()

    @log_function_name
    def _get_vertices_in_components(self) -> list[list[int]]:
        """Reduce the network to the specified component only."""
        vertices_in_components = sorted(nx.connected_components(self.graph), key=len, reverse=True)
        return [list(element) for element in vertices_in_components]

    @log_function_name
    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simplex dimension {simplex_dimension}.'

        num_simplices = self.num_simplices(simplex_dimension)
        degree_sequence: list[int] = []
        for part in self.get_partition(10000):
            degree_sequence_of_part = \
                part.cpp_network.calc_degree_sequence(simplex_dimension, neighbor_dimension)
            degree_sequence.extend(degree_sequence_of_part)

        remaining_simplices = num_simplices - len(degree_sequence)
        degree_sequence.extend([0] * remaining_simplices)

        return sorted(degree_sequence)

    @log_function_name
    def get_partition(
        self,
        min_vertices_in_part: int
    ) -> list[FiniteNetwork]:
        """Return a partition of the network based on the number of vertices."""
        vertices_in_components = self._get_vertices_in_components()
        partitions: list[FiniteNetwork] = []
        interactions = self.interactions  # not to convert interactions any more

        vertices_in_part: list[int] = []
        for vertices_in_component in vertices_in_components:
            vertices_in_part.extend(vertices_in_component)
            if len(vertices_in_part) >= min_vertices_in_part:
                interactions_in_partition = filter_simplices(interactions, vertices_in_part)
                partitions.append(FiniteNetwork(self.max_dimension, vertices_in_part, interactions_in_partition))
                vertices_in_part = []

        # k_cliques = list(nx.community.k_clique_communities(self.graph, max_intersecting_simplex_dimension))
        # vertices_in_part: list[int] = []
        # for k_clique in k_cliques:
        #     vertices_in_part.extend(k_clique)
        #     if len(vertices_in_part) >= min_vertices_in_part:
        #         interactions_in_partition = filter_simplices(interactions, vertices_in_part)
        #         partitions.append(FiniteNetwork(self.max_dimension, vertices_in_part, interactions_in_partition))
        #         vertices_in_part = []

        if len(vertices_in_part) > 0:
            interactions_in_partition = filter_simplices(interactions, vertices_in_part)
            partitions.append(FiniteNetwork(self.max_dimension, vertices_in_part, interactions_in_partition))

        return partitions

    @property
    def components(self) -> list[FiniteNetwork]:
        """Getter of Betti numbers."""
        if self._components is None:
            self._calc_components()
        return self._components

    def _calc_components(self) -> None:
        components: list[FiniteNetwork] = []
        vertices_in_components = self._get_vertices_in_components()
        interactions = self.interactions  # not to convert interactions any more
        for vertices_in_component in vertices_in_components:
            print(f'vertices: {len(vertices_in_component)}', end=' ')
            interactions_in_component = filter_simplices(interactions, vertices_in_component)
            component = FiniteNetwork(self.max_dimension, vertices_in_component, interactions_in_component)
            component.graph = self.graph.subgraph(vertices_in_component)
            components.append(component)
            print('done')
        self._components = components

    @log_function_name
    def _num_edges(self) -> int:
        """Return the number of edges in the graph."""
        return self.graph.number_of_edges()

    def __copy__(self):
        """Shallow copy of self."""
        return self

    def __deepcopy__(self, memo):
        """Deep copy of self."""
        return self
