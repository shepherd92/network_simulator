#!/usr/bin/env python3
"""Generic infinite network type."""

from __future__ import annotations

from copy import deepcopy
from logging import debug, info
from typing import Any

from gudhi.simplex_tree import SimplexTree
import networkx as nx
import numpy as np
import numpy.typing as npt
from tqdm import tqdm

# pylint: disable-next=no-name-in-module
from cpp_modules.build.simplicial_complex import calc_degree_sequence  # type: ignore
from distribution.empirical_distribution import EmpiricalDistribution
from network.network import Network
from network.property import BaseNetworkProperty, DerivedNetworkProperty


class FiniteNetwork(Network):
    """Base class representing a data set or a simulated finite network."""

    def __init__(self, max_dimension: int) -> None:
        """Construct an empty network."""
        super().__init__(max_dimension)
        self._is_persistence_computed: bool = False
        self._is_collapsed_persistence_computed: bool = False
        self._betti_numbers: npt.NDArray[np.int_] | None = None
        self._components: list[FiniteNetwork] | None = None

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        simplicial_complex = SimplexTree()
        for node in tqdm(self.graph.nodes, desc='Inserting nodes from graph to simplical complex', delay=10):
            simplicial_complex.insert((node,))
        for edge in tqdm(self.graph.edges, desc='Inserting edges from graph to simplical complex', delay=10):
            simplicial_complex.insert(edge)
        self.simplicial_complex = simplicial_complex

    def get_component(self, component_index: int) -> FiniteNetwork:
        """Return the network of the specified component."""
        reduced_network = FiniteNetwork(self.max_dimension)
        if component_index != -1:
            components = sorted(nx.connected_components(self.graph), key=len, reverse=True)
            node_ids_in_component = components[component_index]
            reduced_network.graph = self.graph.subgraph(node_ids_in_component)
            reduced_network.simplicial_complex = self.simplicial_complex
            reduced_network.interactions = self.interactions
            reduced_network.filter_to_graph()
            if self.vertex_positions is not None:
                reduced_network.vertex_positions = {
                    node_id: position
                    for node_id, position in self.vertex_positions.items()
                    if node_id in node_ids_in_component
                }
            if self.interaction_positions is not None:
                reduced_network.interaction_positions = {
                    id: position
                    for id, position in self.interaction_positions.items()
                    if id in node_ids_in_component
                }
        return reduced_network

    def reduce_to_component(self, component_index: int) -> None:
        """Reduce the network to the specified component only."""
        if component_index != -1:
            components = sorted(nx.connected_components(self.graph), key=len, reverse=True)
            node_ids_in_component = components[component_index]
            self.graph = self.graph.subgraph(node_ids_in_component)
            self.filter_to_graph()
            if self.vertex_positions is not None:
                self.vertex_positions = {
                    node_id: position
                    for node_id, position in self.vertex_positions.items()
                    if node_id in node_ids_in_component
                }

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

    def collapse(self) -> FiniteNetwork:
        """Collapse edges of the simplicial complex preserving its 1 homology."""
        debug(f'Collapsing network, containing currently: {self.num_vertices} vertices.')
        collapsed_network = FiniteNetwork(self.max_dimension)
        collapsed_network.simplicial_complex = deepcopy(self.simplicial_complex)
        collapsed_network.simplicial_complex.collapse_edges()
        collapsed_network.expand()
        debug(f'Collapsed containing {collapsed_network.num_vertices} vertices.')
        return collapsed_network

    def calc_base_property(self, property_type: BaseNetworkProperty.Type) -> Any:
        """Return a base property of the network."""
        debug(f'Calculating {property_type.name}...')

        if property_type == BaseNetworkProperty.Type.NUM_OF_NODES:
            property_value = self.num_vertices
        elif property_type == BaseNetworkProperty.Type.NUM_OF_EDGES:
            simplex_dimension_value_counts = self.simplex_dimension_distribution.calc_value_counts()
            property_value = simplex_dimension_value_counts[simplex_dimension_value_counts[:, 0] == 1][0, 1] \
                if len(simplex_dimension_value_counts[simplex_dimension_value_counts[:, 0] == 1]) > 0 \
                else 0
        elif property_type == BaseNetworkProperty.Type.NUM_OF_TRIANGLES:
            simplex_dimension_value_counts = self.simplex_dimension_distribution.calc_value_counts()
            if len(simplex_dimension_value_counts[simplex_dimension_value_counts[:, 0] == 2]) == 0:
                property_value = 0
            else:
                property_value = simplex_dimension_value_counts[simplex_dimension_value_counts[:, 0] == 2][0, 1]
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
        elif property_type == BaseNetworkProperty.Type.DIMENSION:
            property_value = self.simplicial_complex.dimension()
        elif property_type == BaseNetworkProperty.Type.NUM_OF_SIMPLICES:
            property_value = self.num_simplices
        elif property_type == BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_interaction_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.SIMPLEX_DIMENSION_DISTRIBUTION:
            property_value = self.simplex_dimension_distribution
        elif property_type == BaseNetworkProperty.Type.FACET_DIMENSION_DISTRIBUTION:
            property_value = self._calculate_facet_dimension_distribution()
        elif property_type == BaseNetworkProperty.Type.BETTI_NUMBERS:
            # property_value = self._calculate_betti_numbers_with_collapse(max_num_of_vertices=1000)
            property_value = self.betti_numbers
        elif property_type == BaseNetworkProperty.Type.BETTI_NUMBERS_BY_COMPONENT:
            property_value = self._calculate_betti_numbers_in_components()
        elif property_type == BaseNetworkProperty.Type.VERTICES_BY_COMPONENT:
            property_value = self._calculate_vertices_in_components()
        elif property_type == BaseNetworkProperty.Type.PERSISTENCE:
            property_value = tuple(self.simplicial_complex.persistence())
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
            'num_of_vertices': self.num_vertices,
            'num_of_interactions': len(self.interactions),
            'num_of_simplices': self.num_simplices,
            'max_dimension': self.max_dimension,
            'num_of_components': len(list(nx.connected_components(self.graph))),
            'num_of_vertices_in_component_0': self.num_of_vertices_in_component(0),
        }

    def _reset(self):
        """Reset saved partial results."""
        super()._reset()
        self._is_persistence_computed = False
        self._betti_numbers = None

    def _calculate_average_degree(self) -> float:
        num_of_nodes = self.graph.number_of_nodes()
        num_of_edges = self.graph.number_of_edges()
        average_degree = num_of_edges / num_of_nodes * 2
        return average_degree

    def _compute_persistence(self) -> None:
        """Compute persistence of the simplicial complex."""
        debug('Computing persistence...')
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

    def _calculate_betti_numbers_in_components(self) -> npt.NDArray[np.int_]:

        betti_numbers_list: list[npt.NDArray[np.int_]] = [
            component.calc_base_property(BaseNetworkProperty.Type.BETTI_NUMBERS)[:, 1]
            for component in tqdm(
                self.components,
                desc='Calculating Betti numbers in components',
                delay=10,
            )
        ]

        betti_numbers_by_component = np.vstack(betti_numbers_list)
        return betti_numbers_by_component

    def _calculate_vertices_in_components(self) -> npt.NDArray[np.int_]:

        vertices_in_components_list: list[int] = [
            component.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)
            for component in tqdm(
                self.components,
                desc='Calculating vertices in components',
                delay=10,
            )
        ]

        betti_numbers_by_component = np.array(vertices_in_components_list)
        return betti_numbers_by_component

    def _calculate_betti_numbers_with_collapse(self, max_num_of_vertices: int) -> npt.NDArray[np.int_]:
        """Calculate the Betti numbers for different dimensions."""
        if self._betti_numbers is None:
            return self._betti_numbers

        collapsed_network = self

        while collapsed_network.num_vertices > max_num_of_vertices:
            debug(f'Collapsing network as it contains {collapsed_network.num_vertices} vertices.')
            collapsed_network = collapsed_network.collapse()

        return collapsed_network._calculate_betti_numbers()  # pylint: disable=protected-access

    def _calculate_betti_numbers(self) -> npt.NDArray[np.int_]:
        """Calculate the Betti numbers for different dimensions."""
        if not self.is_persistence_computed:
            self._compute_persistence()

        if self.simplicial_complex.dimension() != 0:
            betti_numbers = self.simplicial_complex.betti_numbers()
        else:
            betti_numbers = [self.graph.number_of_nodes()]

        betti_number_array = np.zeros((self.max_dimension, 2), dtype=int)
        betti_number_array[:, 0] = range(self.max_dimension)

        for dimension, betti_number in enumerate(betti_numbers):
            betti_number_array[dimension, 1] = betti_number

        return betti_number_array

    def _calc_persistence_pairs(self) -> list[tuple[list[int], list[int]]]:
        if not self.is_persistence_computed:
            self._compute_persistence()
        persistence_pairs = self.simplicial_complex.persistence_pairs()
        return persistence_pairs

    def _get_all_components(self) -> list[FiniteNetwork]:
        """Reduce the network to the specified component only."""
        graph_components = sorted(nx.connected_components(self.graph), key=len, reverse=True)

        components: list[FiniteNetwork] = []
        for graph_component in graph_components:
            component = FiniteNetwork(self.max_dimension)
            component.graph = self.graph.subgraph(graph_component)
            component.generate_clique_complex_from_graph()
            components.append(component)

        return components

    @ property
    def simplices(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        simplices_with_filtration = self.simplicial_complex.get_simplices()
        simplices = [simplex for simplex, _ in simplices_with_filtration]
        return simplices

    @property
    def simplicial_complex(self) -> SimplexTree:
        """Getter of simplicial complex."""
        return super(FiniteNetwork, self).simplicial_complex

    @property
    def betti_numbers(self) -> npt.NDArray[np.int_]:
        """Getter of Betti numbers."""
        if self._betti_numbers is None:
            self._betti_numbers = self._calculate_betti_numbers()
        return self._betti_numbers

    @property
    def components(self) -> list[FiniteNetwork]:
        """Getter of Betti numbers."""
        if self._components is None:
            self._components = self._get_all_components()
        return self._components

    @simplicial_complex.setter
    def simplicial_complex(self, value: SimplexTree) -> None:
        """Setter of simplicial complex."""
        super(FiniteNetwork, self.__class__).simplicial_complex.fset(self, value)
        self._is_persistence_computed = False

    @ property
    def is_persistence_computed(self) -> bool:
        """Return if persistence is computed on the simplicial complex."""
        return self._is_persistence_computed

    def __copy__(self):
        """Shallow copy of self."""
        return self

    def __deepcopy__(self, memo):
        """Deep copy of self."""
        return self

    def __str__(self) -> str:
        """Return a string representation based on the network properties."""
        return '\n'.join([
            f'Number of vertices: {self.num_vertices}',
            f'Number of interactions: {len(self._interactions)}',
            f'Number of simplices: {self.num_simplices}',
            f'Max Dimension: {self.max_dimension}',
            f'Number of components: {len(list(nx.connected_components(self.graph)))}',
            f'Number of vertices in component 0: {self.num_of_vertices_in_component(0)}',
        ])
