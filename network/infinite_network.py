#!/usr/bin/env python3
"""Generic infinite network type."""

from __future__ import annotations

from logging import info
from pathlib import Path
from typing import Any

import pandas as pd

# pylint: disable=no-name-in-module
from cpp_modules.build.cpp_plugin import InfiniteNetwork as CppInfiniteNetwork
# pylint: enable=no-name-in-module
from distribution.empirical_distribution import EmpiricalDistribution
from network.network import Network
from network.property import BaseNetworkProperty
from tools.logging_helper import log_function_name


class InfiniteNetworkSet:
    """A set of infinite network."""

    def __init__(self, infinite_networks: list[InfiniteNetwork]) -> None:
        """Construct an empty network."""
        self._infinite_networks = infinite_networks

    @log_function_name
    def get_largest_network(self) -> InfiniteNetwork | None:
        """Return the largest infinite network."""
        if len(self._infinite_networks) == 0:
            return None

        largest_network = max(
            self._infinite_networks,
            key=lambda network: network.num_simplices(0)
        )
        return largest_network

    @log_function_name
    def calc_network_summary(
        self,
        properties_to_calculate: list[BaseNetworkProperty.Type]
    ) -> dict[BaseNetworkProperty.Type, EmpiricalDistribution]:
        """Calculate the summary of the network."""
        info('Infinite network summary calculation started.')
        summary: dict[BaseNetworkProperty.Type, EmpiricalDistribution] = {
            property_type: self.calc_typical_property_distribution(property_type)
            for property_type in properties_to_calculate
        }
        info('Infinite network summary calculation finished.')
        return summary

    @log_function_name
    def calc_typical_property_distribution(
        self,
        property_type: BaseNetworkProperty.Type,
    ) -> EmpiricalDistribution:
        """Generate typical properties of the given type."""
        typical_property_sets = [
            infinite_network.calc_base_property_value_set(property_type)
            for infinite_network in self._infinite_networks
        ]

        distribution = EmpiricalDistribution([
            value for property_set in typical_property_sets for value in property_set
        ])
        return distribution

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        information = self.get_info_as_dict()
        data_frame = pd.DataFrame(information, index=[0])
        data_frame.to_csv(save_path, index=False)

    def get_info_as_dict(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        return {
            'num_of_networks': len(self._infinite_networks),
            'num_of_simplices': sum([network.num_simplices() for network in self._infinite_networks]),
            'max_dimension': self._infinite_networks[0].max_dimension,
        }


class InfiniteNetwork(Network):
    """Represent an "infinite network" in which network size effects do not play a role."""

    def __init__(self, max_dimension: int, vertices: list[int], interactions: list[list[int]]) -> None:
        """Construct an empty network."""
        super().__init__()
        assert isinstance(max_dimension, int)
        self._cpp_network = CppInfiniteNetwork(max_dimension, vertices, interactions, 0)

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

    @log_function_name
    def _num_edges(self) -> int:
        """Return the number of edges in the graph."""
        return self.graph.degree(0)

    def get_info_as_dict(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        raise NotImplementedError('This method is not implemented.')

    def _calc_typical_in_degree(self) -> list[int]:
        return [self.digraph.in_degree(0)]

    def _calc_typical_out_degree(self) -> list[int]:
        return [self.digraph.out_degree(0)]
