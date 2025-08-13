#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from logging import info

import matplotlib.pyplot as plt

from network.finite_hypergraph import FiniteHypergraph
from network.property import BaseNetworkProperty
from reports.network_analysis.finite_network_analyzer import FiniteNetworkAnalyzer
from reports.plotting_helper import plot_network


class FiniteHypergraphAnalyzer(FiniteNetworkAnalyzer):
    """Class to analyze finite hypergraph models."""

    def analyze_finite_hypergraph(
        self,
        network: FiniteHypergraph,
        calculated_properties: list[BaseNetworkProperty],
    ) -> None:
        """Generate exploratory analysis of the provided finite network."""
        info('Finite network analysis started.')

        plt.rcParams["text.usetex"] = False

        if self._parameters.plot:
            if self._parameters.plot_entire_network:
                plot_network(network, False, self.save_directory / 'network.png')
            if self._parameters.plot_network_giant_component:
                plot_network(network.get_component(0), False, self.save_directory / 'network_giant_component.png')
            if self._parameters.plot_network_determined_positions:
                plot_network(network, True, self.save_directory / 'network_fixed_vertex_positions.png')

        summary = network.calc_network_summary(calculated_properties)

        axes_grid_height = 3
        axes_grid_width = 3

        network.save_info(self.save_directory / 'network_info.csv')
        figure = plt.figure('Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
        axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
        subfigure_row_index = 0

        self.report_vertex_edge_degree_distribution(
            summary.get(BaseNetworkProperty.vertex_edge_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 0])
        )
        self.report_edge_triangle_degree_distribution(
            summary.get(BaseNetworkProperty.edge_triangle_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 1])
        )
        self.report_triangle_tetrahedron_degree_distribution(
            summary.get(BaseNetworkProperty.triangle_tetrahedra_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 2])
        )

        subfigure_row_index += 1

        self.report_simplex_dimension_distribution(
            summary.get(BaseNetworkProperty.simplex_dimension_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 0])
        )

        subfigure_row_index += 1

        self.report_vertices_by_component(
            summary.get(BaseNetworkProperty.num_of_vertices_by_component),
            figure.add_subplot(axes_grid[subfigure_row_index, 0])
        )
        self.report_betti_numbers(
            summary.get(BaseNetworkProperty.betti_numbers),
            figure.add_subplot(axes_grid[subfigure_row_index, 1])
        )
        self.report_betti_number_1_by_component(
            summary.get(BaseNetworkProperty.betti_numbers_by_component),
            figure.add_subplot(axes_grid[subfigure_row_index, 2])
        )

        subfigure_row_index += 1

        figure.tight_layout()
        figure.savefig(self.save_directory / 'whole_report_finite_network.png')
        figure.clf()

        info('Finite network analysis finished.')
