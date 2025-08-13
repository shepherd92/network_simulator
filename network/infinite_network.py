#!/usr/bin/env python3
"""Infinite clique complex network type."""

from logging import info
from pathlib import Path
from typing import Any, Self

import networkx as nx
import pandas as pd

# pylint: disable-next=no-name-in-module
from cpp_plugin.build.release.cpp_plugin import (
    InfiniteCliqueComplex as CppInfiniteCliqueComplex,
    InfiniteHypergraph as CppInfiniteHypergraph
)
from distribution.empirical_distribution import EmpiricalDistribution
from network.network import Network
from network.property import BaseNetworkProperty
from tools.logging_helper import log_function_name


class InfiniteNetwork(Network):
    """Represent an "infinite network" in which network size effects do not play a role."""

    def calc_base_property(self, property_type: BaseNetworkProperty) -> Any:
        """Return a base property of the network.

        Calculate a set of values of a property type.
        """
        if property_type == BaseNetworkProperty.num_of_edges:
            property_value = self.num_simplices(1)
        elif property_type == BaseNetworkProperty.num_of_triangles:
            property_value = self.num_simplices(2)
        elif property_type == BaseNetworkProperty.vertex_edge_degree_distribution:
            property_value = self._calc_degree_sequence(0, 1)
        elif property_type == BaseNetworkProperty.edge_triangle_degree_distribution:
            property_value = self._calc_degree_sequence(1, 2)
        elif property_type == BaseNetworkProperty.triangle_tetrahedra_degree_distribution:
            property_value = self._calc_degree_sequence(2, 3)
        else:
            raise NotImplementedError(
                f'Requested property type {property_type} is not available.'
            )

        return property_value

    def num_all_vertices(self) -> int:
        """Return the number of all vertices."""
        return len(self.cpp_network.get_skeleton(1))

    @log_function_name
    def _generate_graph_from_simplicial_complex(self) -> nx.Graph:
        """Set graph to represent the simplicial complex."""
        # assert self.simplicial_complex.num_vertices() != 0, 'Simplicial complex is not built.'
        skeleton = self.cpp_network.get_skeleton(2)
        graph = nx.Graph()
        graph.add_nodes_from([simplex[0] for simplex in skeleton if len(simplex) == 1])
        graph.add_edges_from([[*simplex] for simplex in skeleton if len(simplex) == 2])
        return graph


class InfiniteNetworkSet:
    """A set of infinite network."""

    def __init__(self, infinite_networks: list[InfiniteNetwork]) -> None:
        """Construct an empty network."""
        self._infinite_networks = infinite_networks

    @log_function_name
    def get_largest_network(self) -> Self | None:
        """Return the largest infinite network."""
        if len(self._infinite_networks) == 0:
            return None

        largest_network = max(
            self._infinite_networks,
            key=lambda network: network.num_all_vertices()
        )
        return largest_network

    @log_function_name
    def calc_network_summary(
        self,
        properties_to_calculate: list[BaseNetworkProperty]
    ) -> dict[BaseNetworkProperty, EmpiricalDistribution]:
        """Calculate the summary of the network."""
        info('Infinite network summary calculation started.')
        summary: dict[BaseNetworkProperty, list[Any]] = {
            property_type: self.calc_base_property(property_type)
            for property_type in properties_to_calculate
        }
        info('Infinite network summary calculation finished.')
        return summary

    @log_function_name
    def calc_base_property(
        self,
        property_type: BaseNetworkProperty,
    ) -> list[Any]:
        """Generate typical properties of the given type."""
        base_property_set = [
            infinite_network.calc_base_property(property_type)
            for infinite_network in self._infinite_networks
        ]
        return base_property_set

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        information = self.info()
        data_frame = pd.DataFrame(information, index=[0])
        data_frame.to_csv(save_path, index=False)

    def info(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        return {
            'num_of_networks': len(self._infinite_networks),
            'max_dimension': self._infinite_networks[0].max_dimension,
        }
