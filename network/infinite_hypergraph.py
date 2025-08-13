#!/usr/bin/env python3
"""Infinite hypergraph network type."""

from typing import Any, NewType

from network.infinite_network import InfiniteNetwork
from network.property import BaseNetworkProperty


CppInfiniteHypergraph = NewType('CppInfiniteHypergraph', Any)


class InfiniteHypergraph(InfiniteNetwork):
    """An "infinite hypergraph network" in which network size effects do not play a role."""

    def __init__(self, cpp_network: CppInfiniteHypergraph) -> None:
        """Construct an empty network."""
        super().__init__(cpp_network)
        self._interaction_positions: dict[int, tuple[float, ...]] | None = None

    def calc_base_property(self, property_type: BaseNetworkProperty) -> Any:
        """Return a base property of the network.

        Calculate a set of values of a property type.
        """
        if property_type == BaseNetworkProperty.vertex_interaction_degree_distribution:
            property_value = self._calc_typical_vertex_interaction_degree()
        elif property_type == BaseNetworkProperty.edge_interaction_degree_distribution:
            property_value = self._calc_typical_edge_interaction_degree()
        elif property_type == BaseNetworkProperty.triangle_interaction_degree_distribution:
            property_value = self._calc_typical_triangle_interaction_degree()
        elif property_type == BaseNetworkProperty.interaction_vertex_degree_distribution:
            property_value = self._calculate_interaction_dimension_distribution()
        else:
            property_value = super().calc_base_property(property_type)

        return property_value

    def _calculate_interaction_dimension_distribution(self) -> list[int]:
        """Return the number of interactions for each dimension."""
        return self.cpp_network.calc_interaction_dimension_distribution()

    def _calc_typical_vertex_interaction_degree(self) -> list[int]:
        """Return the number of interactions for each vertex."""
        return self.cpp_network.calc_simplex_interaction_degree_sequence(0)

    def _calc_typical_edge_interaction_degree(self) -> list[int]:
        """Return the number of interactions for each vertex."""
        return self.cpp_network.calc_simplex_interaction_degree_sequence(1)

    def _calc_typical_triangle_interaction_degree(self) -> list[int]:
        """Return the number of interactions for each vertex."""
        return self.cpp_network.calc_simplex_interaction_degree_sequence(2)

    @property
    def interactions(self) -> list[list[int]]:
        """Get network interactions."""
        return self._cpp_network.interactions

    @interactions.setter
    def interactions(self, value: list[list[int]]) -> None:
        """Setter for interactions."""
        self._cpp_network.interactions = value

    @property
    def interaction_positions(self) -> dict[int, tuple[float, ...]]:
        """Return interaction positions to plot."""
        return self._interaction_positions

    @interaction_positions.setter
    def interaction_positions(self, value: dict[int, tuple[float, ...]]) -> None:
        """Setter of interaction positions."""
        self._interaction_positions = value
