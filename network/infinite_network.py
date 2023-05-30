#!/usr/bin/env python3
"""Generic infinite network type."""

from __future__ import annotations

from gudhi.simplex_tree import SimplexTree
from math import comb

from network.network import Network
from network.property import BaseNetworkProperty


class InfiniteNetwork(Network):
    """Represent an "infinite network" in which network size effects do not play a role."""

    def generate_simplicial_complex_from_graph(self) -> None:
        """Set the simplicial complex to represent the graph."""
        simplicial_complex = SimplexTree()
        self._interactions: list[list[int]] = []
        self._facets: list[list[int]] = []

        for node in self.graph.nodes:
            simplicial_complex.insert((node,))
        for edge in self.graph.edges:
            simplicial_complex.insert(edge)
            if 0 in edge:
                self._interactions.append(edge)

        self.simplicial_complex = simplicial_complex

    def add_simplex(self, simplex: list[int]) -> None:
        """Insert a simplex to the simplicial complex.

        Add the skeleton of the simplex as its dimension is too high.
        """
        if len(simplex) == 0:
            return

        skeleton = self._get_simplex_skeleton_for_max_dimension(simplex)
        self.add_simplices_batch(skeleton)

    def calc_base_property_value_set(self, property_type: BaseNetworkProperty.Type) -> list[float | int]:
        """Return a base property of the network.

        Calculate a set of values of a property type.
        """
        if property_type == BaseNetworkProperty.Type.DEGREE_DISTRIBUTION:
            property_value_set = self._calc_degree_sequence(0, 1)
        elif property_type == BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION:
            property_value_set = self._calc_typical_in_degree()
        elif property_type == BaseNetworkProperty.Type.OUT_DEGREE_DISTRIBUTION:
            property_value_set = self._calc_typical_out_degree()
        elif property_type == BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1:
            property_value_set = self._calc_degree_sequence(1, 2)
        elif property_type == BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2:
            property_value_set = self._calc_degree_sequence(2, 3)
        elif property_type == BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_3:
            property_value_set = self._calc_degree_sequence(3, 4)
        else:
            raise NotImplementedError(
                f'Requested property type {property_type} is not available.'
            )

        return property_value_set

    def _calc_typical_in_degree(self) -> list[int]:
        return self.digraph.in_degree(0)

    def _calc_typical_out_degree(self) -> list[int]:
        return self.digraph.out_degree(0)

    def _calc_degree_sequence(self, simplex_dimension: int, neighbor_dimension: int) -> list[int]:

        assert neighbor_dimension > simplex_dimension, \
            f'Neighbor dimension {neighbor_dimension} must be greater than simlex dimension {simplex_dimension}.'

        # extract facets with dimension (len(facet) - 1) higher or equal to neighbor dimension
        # assumption: every facet contains vertex 0
        facets = filter(lambda facet: len(facet) - 1 >= neighbor_dimension, self.facets)

        degree_sequence: list[int] = []
        for facet in facets:
            facet_dimension = len(facet) - 1
            num_of_simplices_in_facet_containing_0 = comb(facet_dimension, simplex_dimension)
            num_of_neighbors_in_facet_containing_0 = comb(
                facet_dimension - simplex_dimension,
                neighbor_dimension - simplex_dimension
            )
            degree_sequence += [num_of_neighbors_in_facet_containing_0] * num_of_simplices_in_facet_containing_0

        return degree_sequence

    @property
    def simplices(self) -> list[list[int]]:
        """Get the simplices associated to the network."""
        simplices_with_filtration = self.simplicial_complex.get_cofaces([0], 0)
        simplices = [simplex for simplex, _ in simplices_with_filtration]
        return simplices
