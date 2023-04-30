#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from functools import partial
from logging import debug, info
from pathlib import Path
from typing import Any, Callable

from gudhi.persistence_graphical_tools import (
    plot_persistence_barcode,
    plot_persistence_diagram,
)
import matplotlib.pyplot as plt

from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from network.finite_network import FiniteNetwork
from network.infinite_network import InfiniteNetwork
from network.property import BaseNetworkProperty
from reports.plotting_helper import (
    approximate_and_plot_pdf,
    plot_empirical_distribution_histogram_with_info,
    plot_finite_network,
    plot_value_counts,
    print_not_calculated,
    print_info,
)


def analyze_finite_network(
    network: FiniteNetwork,
    calculated_properties: list[BaseNetworkProperty.Type],
    save_directory: Path
) -> None:
    """Generate exploratory analysis of the provided finite network."""
    info('Finite network analysis started.')

    plt.rcParams["text.usetex"] = False

    axes_grid_height = 7
    axes_grid_width = 3

    figure = plt.figure("Network Analysis", figsize=(axes_grid_width * 10, axes_grid_height * 10))
    axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)

    network_to_analyze = network

    debug('Summary calculation started.')
    summary = network_to_analyze.calc_network_summary(calculated_properties)
    debug('Summary calculation finished.')

    debug('Plotting simplicial complex started.')
    network_to_plot = network_to_analyze.get_component(0)
    simplicial_complex_axes = figure.add_subplot(axes_grid[0:axes_grid_width, 0:axes_grid_width])
    plot_finite_network(network_to_plot, simplicial_complex_axes)
    figure.savefig(save_directory / 'network.png')
    print_info([network], simplicial_complex_axes)
    debug('Plotting simplicial complex finished.')

    subfigure_row_index = axes_grid_width

    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.DEGREE_DISTRIBUTION, None),
        'Total Degree Distribution',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POWER_LAW),
        figure.add_subplot(axes_grid[subfigure_row_index, 0])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION, None),
        'In Degree Distribution',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POWER_LAW),
        figure.add_subplot(axes_grid[subfigure_row_index, 1])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.OUT_DEGREE_DISTRIBUTION, None),
        'Out Degree Distribution',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POISSON),
        figure.add_subplot(axes_grid[subfigure_row_index, 2])
    )
    subfigure_row_index += 1

    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1, None),
        'Higher-Order Degree Distribution - Dimension 1',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POWER_LAW),
        figure.add_subplot(axes_grid[subfigure_row_index, 0])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2, None),
        'Higher-Order Degree Distribution - Dimension 2',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POWER_LAW),
        figure.add_subplot(axes_grid[subfigure_row_index, 1])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_3, None),
        'Higher-Order Degree Distribution - Dimension 3',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POWER_LAW),
        figure.add_subplot(axes_grid[subfigure_row_index, 2])
    )
    subfigure_row_index += 1

    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.SIMPLEX_DIMENSION_DISTRIBUTION, None),
        'Simplex dimension distribution',
        partial(
            plot_empirical_distribution_histogram_with_info,
            histogram_type=EmpiricalDistribution.HistogramType.INTEGERS
        ),
        figure.add_subplot(axes_grid[subfigure_row_index, 0])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.FACET_DIMENSION_DISTRIBUTION, None),
        'Facet dimension distribution',
        partial(
            plot_empirical_distribution_histogram_with_info,
            histogram_type=EmpiricalDistribution.HistogramType.INTEGERS
        ),
        figure.add_subplot(axes_grid[subfigure_row_index, 1])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION, None),
        'Interaction dimension distribution',
        partial(
            plot_empirical_distribution_histogram_with_info,
            histogram_type=EmpiricalDistribution.HistogramType.INTEGERS
        ),
        figure.add_subplot(axes_grid[subfigure_row_index, 2])
    )
    subfigure_row_index += 1

    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.BETTI_NUMBERS, None),
        'Betti Numbers',
        plot_value_counts,
        figure.add_subplot(axes_grid[subfigure_row_index, 0])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.PERSISTENCE, None),
        'Persistence Diagram',
        plot_persistence_diagram,
        figure.add_subplot(axes_grid[subfigure_row_index, 1])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.PERSISTENCE, None),
        'Persistence Barcode',
        plot_persistence_barcode,
        figure.add_subplot(axes_grid[subfigure_row_index, 2])
    )
    subfigure_row_index += 1

    figure.tight_layout()
    figure.savefig(save_directory / 'whole_report.png')
    figure.clf()

    info('Finite network analysis finished.')


def analyze_infinite_network(
    network: InfiniteNetwork,
    calculated_properties: list[BaseNetworkProperty.Type],
    save_directory: Path
) -> None:
    """Analyze the given infinite network."""
    a = 1
    pass


def _generate_descriptive_report(summary: dict[BaseNetworkProperty.Type, Any]) -> None:

    with open('../data_analysis/data_report.txt', 'w', encoding='utf8') as file:
        file.write(
            f'Number of nodes: {summary.num_of_nodes}\n' +
            f'Number of edges: {summary.num_of_edges}\n' +
            f'Maximum degree: {summary.max_degree}\n' +
            f'Average clustering: {summary.avg_clustering}\n' +
            f'Number of connected components: {summary.num_of_connected_components}\n' +
            f'Dimension: {summary.dimension}\n' +
            f'Number of simplices: {summary.num_of_simplices}\n' +
            f'Betti numbers: {summary.betti_numbers}\n'
        )


def _plot_base_property(
    property_value: Any,
    title: str,
    plotter: Callable[[Any, plt.Axes], None],
    axes: plt.Axes
) -> None:

    axes.set_title(title)

    if property_value is None:
        print_not_calculated(axes)
        return

    plotter(property_value, axes=axes)
