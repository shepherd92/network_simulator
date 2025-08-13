#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from dataclasses import dataclass
from logging import info

import matplotlib.pyplot as plt

from network.infinite_hypergraph import InfiniteHypergraph
from network.infinite_network import InfiniteNetworkSet
from network.property import BaseNetworkProperty
from reports.network_analysis.network_analyzer import NetworkAnalyzer
from reports.plotting_helper import plot_network


class InfiniteNetworkAnalyzer(NetworkAnalyzer):
    """Class to analyze finite hypergraph models."""

    @dataclass
    class Parameters(NetworkAnalyzer.Parameters):
        """Parameters for the network analysis."""
        plot_entire_network: bool
        plot_network_determined_positions: bool

    def __init__(self, parameters: Parameters):
        """Initialize the network analyzer with given parameters."""
        self._parameters = parameters

    def plot_infinite_network(self, network: InfiniteHypergraph):

        if self._parameters.plot_entire_network:
            plot_network(network, False, self.save_directory / 'infinite_network.png')
        if self._parameters.plot_network_determined_positions:
            plot_network(network, True, self.save_directory / 'infinite_network_fixed_vertex_positions.png')

    def analyze_infinite_network_set(
        self,
        network: InfiniteHypergraph | InfiniteNetworkSet,
        calculated_properties: list[BaseNetworkProperty]
    ) -> None:
        """Analyze the given infinite network set."""
        raise NotImplementedError
