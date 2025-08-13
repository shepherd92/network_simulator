#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from dataclasses import dataclass
from logging import info
from pathlib import Path

import matplotlib.pyplot as plt
from network.infinite_hypergraph import InfiniteHypergraph
from network.infinite_network import InfiniteNetworkSet
from network.property import BaseNetworkProperty
from reports.network_analyzer import NetworkAnalyzer
from reports.plotting_helper import (
    plot_network,
)


class InfiniteHypergraphAnalyzer(NetworkAnalyzer):
    """Class to analyze finite hypergraph models."""

    @dataclass
    class Parameters(NetworkAnalyzer.Parameters):
        """Parameters for the network analysis."""
        plot_entire_network: bool
        plot_network_determined_positions: bool

    def __init__(self, parameters: Parameters):
        """Initialize the network analyzer with given parameters."""
        self._parameters = parameters

    def analyze_infinite_hypergraph_set(
        self,
        network: InfiniteNetworkSet,
        calculated_properties: list[BaseNetworkProperty]
    ) -> None:
        """Analyze the given infinite network set."""
        info('Infinite network set analysis started.')

        self.analyze_infinite_hypergraph(network, calculated_properties)

        info('Infinite network set analysis finished.')

    def analyze_infinite_hypergraph(
        self,
        network: InfiniteHypergraph | InfiniteNetworkSet,
        calculated_properties: list[BaseNetworkProperty]
    ) -> None:
        """Analyze the given infinite network set."""
        info('Infinite hypergraph analysis started.')

        plt.rcParams["text.usetex"] = False

        if self._parameters.plot:
            if self._parameters.plot_entire_network:
                plot_network(network, False, self.save_directory / 'infinite_network.png')
            if self._parameters.plot_network_determined_positions:
                plot_network(network, True, self.save_directory / 'infinite_network_fixed_vertex_positions.png')

        summary = network.calc_network_summary(calculated_properties)

        axes_grid_height = 4
        axes_grid_width = 3

        network.save_info(self.save_directory / 'network_info.csv')

        figure = plt.figure('Infinite Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
        axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
        subfigure_row_index = 0

        self.report_vertex_edge_degree_distribution(
            summary.get(BaseNetworkProperty.vertex_edge_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        )
        self.report_edge_triangle_degree_distribution(
            summary.get(BaseNetworkProperty.edge_triangle_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        )
        self.report_triangle_tetrahedron_degree_distribution(
            summary.get(BaseNetworkProperty.triangle_tetrahedra_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        )

        subfigure_row_index += 1

        figure.tight_layout()
        figure.savefig(self.save_directory / 'whole_report_infinite_network.png')
        figure.clf()

        info('Infinite hypergraph analysis finished.')
