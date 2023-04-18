#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from functools import partial
from logging import debug, info
from operator import itemgetter
from pathlib import Path
from typing import Any, Callable

from gudhi.persistence_graphical_tools import (
    plot_persistence_barcode,
    plot_persistence_diagram,
)
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.spatial import ConvexHull

from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from network.finite_network import FiniteNetwork
from network.property import BaseNetworkProperty
from reports.plotting_helper import (
    plot_empirical_distribution_histogram_with_info,
    plot_value_counts,
    approximate_and_plot_pdf,
    print_not_calculated,
    print_info
)


def analyze_network(
    network: FiniteNetwork,
    calculated_properties: list[BaseNetworkProperty.Type],
    path: Path
) -> None:
    """Generate exploratory analysis of the provided network."""
    info('Network analysis started.')

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
    _plot_simplicial_complex(network_to_plot, simplicial_complex_axes)
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
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTIONS, {}).get(1, None),
        'Higher-Order Degree Distribution - Dimension 1',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POWER_LAW),
        figure.add_subplot(axes_grid[subfigure_row_index, 0])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTIONS, {}).get(2, None),
        'Higher-Order Degree Distribution - Dimension 2',
        partial(approximate_and_plot_pdf, theoretical_distribution_type=TheoreticalDistribution.Type.POWER_LAW),
        figure.add_subplot(axes_grid[subfigure_row_index, 1])
    )
    _plot_base_property(
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTIONS, {}).get(3, None),
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
    figure.savefig(path)
    figure.clf()

    info('Network analysis finished.')


def _plot_simplicial_complex(network: FiniteNetwork, axes: plt.Axes) -> None:

    def get_face_color(simplex: set[int]):
        simplex_colors = {
            0: 'black',
            1: 'black',
            2: plt.cm.Blues(0.6),
            3: plt.cm.Purples(0.6),
            4: plt.cm.Oranges(0.6),
            5: plt.cm.Reds(0.6),
        }

        max_dimension_can_be_colored = list(simplex_colors)[-1]
        dimension = len(simplex) - 1
        face_color = simplex_colors[dimension] \
            if dimension <= max_dimension_can_be_colored \
            else simplex_colors[max_dimension_can_be_colored]
        return face_color

    if network.graph.number_of_nodes() <= 1000:
        graphviz_method = 'neato'  # best for a few nodes
    else:
        graphviz_method = 'sfdp'  # scalable method
    node_positions = nx.nx_agraph.graphviz_layout(network.graph, prog=graphviz_method)

    all_facets = network.extract_facets()
    # remove facets with dimension 0 and 1 (points and edges)
    facets_to_plot = sorted([facet for facet in all_facets if len(facet) - 1 > 1], key=len)

    debug(f'Number of facets to plot: {len(facets_to_plot)}')

    simplex_node_positions = [
        list(itemgetter(*facet)(node_positions))
        for facet in facets_to_plot
    ]
    convex_hulls = [
        ConvexHull(simplex_node_pos, qhull_options='QJ Pp')
        for simplex_node_pos in simplex_node_positions
    ]

    polygon_coordinates = [
        np.array([
            node_positions[node_index]
            for node_index in convex_hull.vertices
        ])
        for convex_hull, node_positions in zip(convex_hulls, simplex_node_positions)
    ]
    face_colors = list(map(get_face_color, facets_to_plot))

    polygon_collection = PolyCollection(
        polygon_coordinates,
        facecolors=face_colors,
        edgecolors=('black',),
        linewidths=(0.01,),
        alpha=(0.3,),
    )

    axes.add_collection(polygon_collection)

    nx.draw_networkx_edges(network.graph, node_positions, ax=axes, edge_color='black', width=0.01, alpha=0.3)
    nx.draw_networkx_nodes(network.graph, node_positions, ax=axes, node_color='black', node_size=0.1)

    axes.set_title("Simplicial Complex")
    axes.set_axis_off()


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
