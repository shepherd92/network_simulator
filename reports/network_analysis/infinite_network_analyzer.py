#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from dataclasses import dataclass
from logging import info
from typing import Any

import matplotlib.pyplot as plt

from distribution.empirical_distribution import EmpiricalDistribution
from network.infinite_hypergraph import InfiniteHypergraph
from network.infinite_network import InfiniteNetworkSet
from network.property import BaseNetworkProperty
from reports.network_analysis.network_analyzer import NetworkAnalyzer
from reports.plotting_helper import plot_hypergraph


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

    def plot(self, network: InfiniteHypergraph):

        if self._parameters.plot_entire_network:
            plot_hypergraph(network, False, self.save_directory / 'infinite_network.png')
        if self._parameters.plot_network_determined_positions:
            plot_hypergraph(network, True, self.save_directory / 'infinite_network_fixed_vertex_positions.png')

    def analyze_infinite_network_set(
        self,
        network_set: InfiniteNetworkSet,
        calculated_properties: list[BaseNetworkProperty]
    ) -> None:
        """Analyze the given infinite network set."""
        info('Infinite hypergraph analysis started.')

        plt.rcParams["text.usetex"] = False

        summary = network_set.calc_network_summary(calculated_properties)

        axes_grid_height = 2
        axes_grid_width = 3

        network_set.save_info(self.save_directory / 'network_info.csv')

        figure = plt.figure('Infinite Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
        axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
        subfigure_row_index = 0

        distribution = self._flatten_distribution(summary.get(BaseNetworkProperty.vertex_edge_degree_distribution))
        self.report_vertex_edge_degree_distribution(
            distribution,
            figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        )
        distribution = self._flatten_distribution(summary.get(BaseNetworkProperty.edge_triangle_degree_distribution))
        self.report_edge_triangle_degree_distribution(
            distribution,
            figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        )
        distribution = self._flatten_distribution(summary.get(BaseNetworkProperty.triangle_tetrahedra_degree_distribution))
        self.report_triangle_tetrahedron_degree_distribution(
            distribution,
            figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        )

        subfigure_row_index += 1

        distribution = self._flatten_distribution(summary.get(BaseNetworkProperty.vertex_interaction_degree_distribution))
        self.report_vertex_interaction_degree_distribution(
            distribution,
            figure.add_subplot(axes_grid[subfigure_row_index, 0])
        )
        distribution = self._flatten_distribution(summary.get(BaseNetworkProperty.edge_interaction_degree_distribution))
        self.report_edge_interaction_degree_distribution(
            distribution,
            figure.add_subplot(axes_grid[subfigure_row_index, 1])
        )
        distribution = self._flatten_distribution(summary.get(BaseNetworkProperty.triangle_interaction_degree_distribution))
        self.report_interaction_dimension_distribution(
            distribution,
            figure.add_subplot(axes_grid[subfigure_row_index, 2])
        )
    
        figure.tight_layout()
        figure.savefig(self.save_directory / 'whole_report_infinite_network.png')
        figure.clf()

        info('Infinite hypergraph analysis finished.')

    @staticmethod
    def _flatten_distribution(distribution: list[list[Any]] | None) -> EmpiricalDistribution:
        """Flatten a nested list of integers."""
        return EmpiricalDistribution([value for element in distribution for value in element]) \
            if distribution is not None else None