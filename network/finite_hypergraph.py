#!/usr/bin/env python3
"""Finite hypergraph network type."""

from logging import debug
from typing import Any, Self

# pylint: disable-next=no-name-in-module
from cpp_plugin.build.release.cpp_plugin import FiniteHypergraph as CppFiniteHypergraph
from distribution.empirical_distribution import EmpiricalDistribution
from network.finite_network import FiniteNetwork
from network.property import BaseNetworkProperty
from tools.logging_helper import log_function_name


class FiniteHypergraph(FiniteNetwork):
    """Base class representing a data set or a simulated finite network."""

    def __init__(self, cpp_network: CppFiniteHypergraph) -> None:
        """Construct an empty network."""
        super().__init__(cpp_network)
        self._interaction_positions: dict[int, tuple[float, ...]] | None = None

    @log_function_name
    def calc_base_property(self, property_type: BaseNetworkProperty) -> Any:
        """Return a base property of the network."""
        debug(f'Calculating {property_type.name}...')

        if property_type == BaseNetworkProperty.num_of_interactions:
            property_value = len(self.interactions)
        elif property_type == BaseNetworkProperty.mean_degree:
            property_value = self._calculate_average_degree()
        elif property_type == BaseNetworkProperty.vertex_interaction_degree_distribution:
            property_value = self._calculate_interaction_degree_distribution(0)
        elif property_type == BaseNetworkProperty.edge_interaction_degree_distribution:
            property_value = self._calculate_edge_interaction_degree_distribution()
        elif property_type == BaseNetworkProperty.triangle_interaction_degree_distribution:
            property_value = self._calculate_triangle_interaction_degree_distribution()
        elif property_type == BaseNetworkProperty.interaction_vertex_degree_distribution:
            property_value = self._calculate_interaction_dimension_distribution()
            property_value.value_sequence += 1
        elif property_type == BaseNetworkProperty.persistence_intervals:
            property_value = self._calc_persistence_intervals()
        elif property_type == BaseNetworkProperty.persistence_pairs:
            property_value = self._calc_persistence_pairs()
        else:
            property_value = super().calc_base_property(property_type)

        debug('finished.')

        return property_value

    def info(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        return super().info() | {
            'num_of_interactions': len(self.interactions),
        }

    @log_function_name
    def _calculate_vertex_interaction_degree_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each vertex."""
        return EmpiricalDistribution(self.cpp_network.calc_simplex_interaction_degree_sequence(0))

    @log_function_name
    def _calculate_edge_interaction_degree_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each vertex."""
        return EmpiricalDistribution(self.cpp_network.calc_simplex_interaction_degree_sequence(1))

    @log_function_name
    def _calculate_triangle_interaction_degree_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each vertex."""
        return EmpiricalDistribution(self.cpp_network.calc_simplex_interaction_degree_sequence(2))

    @log_function_name
    def _calculate_interaction_dimension_distribution(self) -> EmpiricalDistribution:
        """Return the number of interactions for each dimension."""
        return EmpiricalDistribution(self.cpp_network.calc_interaction_dimension_distribution())

    @log_function_name
    def _calculate_interaction_degree_distribution(self, simplex_dimension: int) -> EmpiricalDistribution:

        interaction_degree_sequence: list[int] = []
        for part in self.get_partition(10000):
            interaction_degree_sequence_of_part = \
                part.cpp_network.calc_simplex_interaction_degree_sequence(simplex_dimension)
            interaction_degree_sequence.extend(interaction_degree_sequence_of_part)

        return EmpiricalDistribution(interaction_degree_sequence)

    @log_function_name
    def _calc_persistence_intervals(self) -> list[tuple[list[int], list[int]]]:
        return self.cpp_network.calc_persistence_intervals()

    @log_function_name
    def _calc_persistence_pairs(self) -> list[tuple[list[int], list[int]]]:
        return self.cpp_network.calc_persistence_pairs()

    def _construct_self_from_cpp(self, cpp_network: CppFiniteHypergraph) -> Self:
        """Construct self from the given cpp network."""
        return FiniteHypergraph(cpp_network)

    @property
    def weighted(self) -> bool:
        """Return if the cpp_network is weighted."""
        return self.cpp_network.is_weighted()

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
