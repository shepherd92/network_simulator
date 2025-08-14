#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from dataclasses import dataclass
from logging import info

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd

from network.finite_clique_complex import FiniteCliqueComplex
from network.finite_hypergraph import FiniteHypergraph
from network.property import BaseNetworkProperty
from reports.network_analysis.network_analyzer import NetworkAnalyzer
from reports.plotting_helper import (
    check_calculated,
    plot_hypergraph,
    plot_persistence_barcode_,
    plot_value_counts,
    plot_value_counts_log,
)


class FiniteNetworkAnalyzer(NetworkAnalyzer):
    """Class to analyze finite hypergraph models."""

    @dataclass
    class Parameters(NetworkAnalyzer.Parameters):
        """Parameters for the network analysis."""
        plot_entire_network: bool
        plot_network_giant_component: bool
        plot_network_determined_positions: bool

    def __init__(self, parameters: Parameters):
        """Initialize the network analyzer with given parameters."""
        self._parameters = parameters

    def plot(self, network: FiniteHypergraph):

        if self._parameters.plot_entire_network:
            plot_hypergraph(network, False, self.save_directory / 'network.png')
        if self._parameters.plot_network_giant_component:
            plot_hypergraph(network.get_component(0), False, self.save_directory / 'network_giant_component.png')
        if self._parameters.plot_network_determined_positions:
            plot_hypergraph(network, True, self.save_directory / 'network_fixed_vertex_positions.png')

    def analyze_finite_hypergraph(
        self,
        network: FiniteHypergraph | FiniteCliqueComplex,
        calculated_properties: list[BaseNetworkProperty],
    ) -> None:
        """Generate exploratory analysis of the provided finite network."""
        info('Finite network analysis started.')

        plt.rcParams["text.usetex"] = False

        summary = network.calc_network_summary(calculated_properties)

        axes_grid_height = 4
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

        self.report_vertex_interaction_degree_distribution(
            summary.get(BaseNetworkProperty.vertex_interaction_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 0])
        )
        self.report_edge_interaction_degree_distribution(
            summary.get(BaseNetworkProperty.edge_interaction_degree_distribution),
            figure.add_subplot(axes_grid[subfigure_row_index, 1])
        )
        self.report_interaction_dimension_distribution(
            summary.get(BaseNetworkProperty.interaction_vertex_degree_distribution),
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
        self.report_persistence_diagram(
            summary.get(BaseNetworkProperty.persistence_intervals),
            figure.add_subplot(axes_grid[subfigure_row_index, 2])
        )
        # self.report_betti_number_1_by_component(
        #     summary.get(BaseNetworkProperty.betti_numbers_by_component),
        #     figure.add_subplot(axes_grid[subfigure_row_index, 2])
        # )

        subfigure_row_index += 1

        figure.tight_layout()
        figure.savefig(self.save_directory / 'whole_report_finite_network.png')
        figure.clf()

        info('Finite network analysis finished.')

    @check_calculated('Betti numbers')
    def report_betti_numbers(
        self,
        betti_numbers: list[int],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report Betti numbers."""
        betti_numbers_array = np.c_[np.arange(len(betti_numbers)), np.array(betti_numbers)]
        plot_value_counts(betti_numbers_array, axes)
        pd.DataFrame(
            betti_numbers_array,
            columns=['dimension', 'betti_number'],
            dtype=np.int32,
        ).to_csv(self.save_directory / 'betti_numbers.csv', index=False)

    @check_calculated('Betti number 1 by component')
    def report_betti_number_1_by_component(
        self,
        betti_numbers_by_component: list[list[int]],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report Betti number 1 by component."""
        betti_numbers_1 = np.array([
            betti_numbers_in_component[1]
            for betti_numbers_in_component in betti_numbers_by_component
        ])

        values_to_plot = np.c_[np.arange(len(betti_numbers_1)), betti_numbers_1,]
        plot_value_counts_log(values_to_plot, axes)
        data_frame = pd.DataFrame(
            betti_numbers_by_component,
            columns=list(range(len(betti_numbers_by_component[0]))),
            dtype=np.int32,
        )
        data_frame.index.name = 'component_index'
        data_frame.to_csv(self.save_directory / 'betti_numbers_by_component.csv')

    @check_calculated('Vertices by component')
    def report_vertices_by_component(
        self,
        vertices_by_component: npt.NDArray[np.int_],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report vertices by component."""
        values_to_plot = np.c_[
            np.arange(1, len(vertices_by_component) + 1),
            vertices_by_component,
        ]
        plot_value_counts_log(values_to_plot, axes)
        pd.DataFrame(
            values_to_plot,
            columns=['component_index', 'num_of_vertices'],
            dtype=np.int32,
        ).to_csv(self.save_directory / 'vertices_by_component.csv', index=False)

    @check_calculated('Persistence diagram')
    def report_persistence_diagram(
        self,
        persistence: list[list[tuple[float, float]]],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report persistence diagram."""
        self.plot_persistence_diagram_(persistence, axes)
        self._save_persistence(persistence)

    @check_calculated('Persistence barcode')
    def report_persistence_barcode(
        self,
        persistence: list[list[tuple[float, float]]],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report persistence barcode."""
        plot_persistence_barcode_(persistence, axes)
        self._save_persistence(persistence)

    def _save_persistence(self, persistence: list[list[tuple[float, float]]]) -> None:
        persistence_df = pd.DataFrame([
                [dimension, birth, death]
                for dimension, intervals in enumerate(persistence)
                for birth, death in intervals
            ],
            columns=['dimension', 'birth', 'death']
        ).astype(int)

        persistence_df = persistence_df.groupby(persistence_df.columns.tolist(), as_index=False).size()

        persistence_df.to_csv(self.save_directory / 'persistence.csv', index=False)
