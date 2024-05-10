#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from logging import info
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from config_files.distribution_fitting_params import POWER_LAW_FITTING_MINIMUM_VALUE_MODEL

from reports.finite_network_analysis import (
    report_betti_number_1_by_component,
    report_betti_numbers,
    report_edge_triangle_degree_distribution,
    report_facet_dimension_distribution,
    report_in_degree_distribution,
    report_interaction_dimension_distribution,
    report_out_degree_distribution,
    report_simplex_dimension_distribution,
    report_triangle_tetrahedron_degree_distribution,
    report_vertex_edge_degree_distribution,
    report_vertex_interaction_degree_distribution,
    report_vertices_by_component,
    report_persistence_diagram,
)
from network.finite_network import FiniteNetwork
from network.infinite_network import InfiniteNetworkSet
from network.property import BaseNetworkProperty
from reports.plotting_helper import (
    plot_giant_component,
    plot_hypergraph_determined_positions,
    plot_network_determined_positions,
)


PLOT_HYPERGRAPH_DETERMINED_POSITIONS = False


def analyze_model_example_finite_network(
    network: FiniteNetwork,
    calculated_properties: list[BaseNetworkProperty],
    plot: bool,
    save_directory: Path,
) -> None:
    """Generate exploratory analysis of the provided finite network."""
    info('Finite network analysis started.')

    plt.rcParams["text.usetex"] = False

    if plot:
        plot_giant_component(network, save_directory / 'network.png')
        plot_network_determined_positions(network, save_directory / 'network_fixed_vertex_positions.png')
        if PLOT_HYPERGRAPH_DETERMINED_POSITIONS:
            plot_hypergraph_determined_positions(network, save_directory / 'hypergraph_fixed_positions.png')

    summary = network.calc_network_summary(calculated_properties)

    axes_grid_height = 4
    axes_grid_width = 3

    network.save_info(save_directory / 'network_info.csv')
    if BaseNetworkProperty.edges in calculated_properties:
        edges = summary.get(BaseNetworkProperty.edges)
        np.savetxt(save_directory / 'edges.csv', edges, delimiter=',', fmt='%i')

    figure = plt.figure('Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
    axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
    subfigure_row_index = 0

    report_vertex_edge_degree_distribution(
        summary.get(BaseNetworkProperty.vertex_edge_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_in_degree_distribution(
        summary.get(BaseNetworkProperty.in_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_out_degree_distribution(
        summary.get(BaseNetworkProperty.out_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )

    subfigure_row_index += 1

    report_vertices_by_component(
        summary.get(BaseNetworkProperty.num_of_vertices_by_component),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
    )
    report_edge_triangle_degree_distribution(
        summary.get(BaseNetworkProperty.edge_triangle_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_triangle_tetrahedron_degree_distribution(
        summary.get(BaseNetworkProperty.triangle_tetrahedra_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )

    subfigure_row_index += 1

    report_simplex_dimension_distribution(
        summary.get(BaseNetworkProperty.simplex_dimension_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_facet_dimension_distribution(
        summary.get(BaseNetworkProperty.facet_dimension_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_vertex_interaction_degree_distribution(
        summary.get(BaseNetworkProperty.vertex_interaction_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_interaction_dimension_distribution(
        summary.get(BaseNetworkProperty.interaction_vertex_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )

    subfigure_row_index += 1

    report_betti_numbers(
        summary.get(BaseNetworkProperty.betti_numbers),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory
    )
    report_persistence_diagram(
        summary.get(BaseNetworkProperty.persistence_intervals),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory
    )
    report_betti_number_1_by_component(
        summary.get(BaseNetworkProperty.betti_numbers_by_component),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )

    subfigure_row_index += 1

    figure.tight_layout()
    figure.savefig(save_directory / 'whole_report_finite_network.png')
    figure.clf()

    info('Finite network analysis finished.')


def analyze_model_example_infinite_network_set(
    network_set: InfiniteNetworkSet,
    calculated_properties: list[BaseNetworkProperty],
    save_directory: Path
) -> None:
    """Analyze the given infinite network set."""
    info('Infinite network set analysis started.')

    plt.rcParams["text.usetex"] = False

    summary = network_set.calc_network_summary(calculated_properties)
    network_with_most_nodes = network_set.get_largest_network()
    assert network_with_most_nodes is not None, \
        'Infinite network set to be analyzed has no networks.'
    plot_giant_component(network_with_most_nodes, save_directory / 'largest_infinite_network.png')

    axes_grid_height = 4
    axes_grid_width = 3

    network_set.save_info(save_directory / 'network_info.csv')

    figure = plt.figure('Infinite Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
    axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
    subfigure_row_index = 0

    report_vertex_edge_degree_distribution(
        summary.get(BaseNetworkProperty.vertex_edge_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_in_degree_distribution(
        summary.get(BaseNetworkProperty.in_degree_distribution),
        POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory
    )
    report_out_degree_distribution(
        summary.get(BaseNetworkProperty.out_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )

    subfigure_row_index += 1

    report_edge_triangle_degree_distribution(
        summary.get(BaseNetworkProperty.edge_triangle_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_triangle_tetrahedron_degree_distribution(
        summary.get(BaseNetworkProperty.triangle_tetrahedra_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )

    subfigure_row_index += 1

    figure.tight_layout()
    figure.savefig(save_directory / 'whole_report_infinite_networks.png')
    figure.clf()

    info('Infinite network analysis finished.')
