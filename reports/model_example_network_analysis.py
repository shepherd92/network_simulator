#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from logging import info
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from config_files.distribution_fitting_params import POWER_LAW_FITTING_MINIMUM_VALUE_MODEL

from reports.finite_network_analysis import (
    # report_betti_number_1_by_component,
    report_betti_numbers,
    report_edge_interaction_degree_distribution,
    report_edge_triangle_degree_distribution,
    report_interaction_dimension_distribution,
    report_simplex_dimension_distribution,
    report_triangle_tetrahedron_degree_distribution,
    report_vertex_edge_degree_distribution,
    report_vertex_interaction_degree_distribution,
    report_vertices_by_component,
    report_persistence_diagram,
)
from network.finite_hypergraph import FiniteHypergraph
from network.infinite_hypergraph import InfiniteHypergraph, InfiniteHypergraphSet
from network.property import BaseNetworkProperty
from reports.plotting_helper import (
    plot_hypergraph,
    plot_network,
)


PLOT_ENTIRE_NETWORK = False
PLOT_NETWORK_GIANT_COMPONENT = True
PLOT_NETWORK_DETERMINED_POSITIONS = True
PLOT_HYPERGRAPH_DETERMINED_POSITIONS = False


def analyze_model_example_finite_network(
    network: FiniteHypergraph,
    calculated_properties: list[BaseNetworkProperty],
    plot: bool,
    save_directory: Path,
) -> None:
    """Generate exploratory analysis of the provided finite network."""
    info('Finite network analysis started.')

    plt.rcParams["text.usetex"] = False

    if plot:
        if PLOT_ENTIRE_NETWORK:
            plot_network(network, False, save_directory / 'network.png')
        if PLOT_NETWORK_GIANT_COMPONENT:
            plot_network(network.get_component(0), False, save_directory / 'network_giant_component.png')
        if PLOT_NETWORK_DETERMINED_POSITIONS:
            plot_network(network, True, save_directory / 'network_fixed_vertex_positions.png')
        if PLOT_HYPERGRAPH_DETERMINED_POSITIONS:
            plot_hypergraph(network, True, save_directory / 'hypergraph_fixed_positions.png')

    summary = network.calc_network_summary(calculated_properties)

    axes_grid_height = 5
    axes_grid_width = 3

    network.save_info(save_directory / 'network_info.csv')
    figure = plt.figure('Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
    axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
    subfigure_row_index = 0

    subfigure_row_index += 1

    report_vertex_edge_degree_distribution(
        summary.get(BaseNetworkProperty.vertex_edge_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
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

    report_vertex_interaction_degree_distribution(
        summary.get(BaseNetworkProperty.vertex_interaction_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )
    report_edge_interaction_degree_distribution(
        summary.get(BaseNetworkProperty.edge_interaction_degree_distribution),
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

    report_simplex_dimension_distribution(
        summary.get(BaseNetworkProperty.simplex_dimension_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
    )

    subfigure_row_index += 1

    report_vertices_by_component(
        summary.get(BaseNetworkProperty.num_of_vertices_by_component),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
    )
    report_betti_numbers(
        summary.get(BaseNetworkProperty.betti_numbers),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory
    )
    report_persistence_diagram(
        summary.get(BaseNetworkProperty.persistence_intervals),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )
    # report_betti_number_1_by_component(
    #     summary.get(BaseNetworkProperty.betti_numbers_by_component),
    #     figure.add_subplot(axes_grid[subfigure_row_index, 2]),
    #     save_directory
    # )

    subfigure_row_index += 1

    figure.tight_layout()
    figure.savefig(save_directory / 'whole_report_finite_network.png')
    figure.clf()

    info('Finite network analysis finished.')


def analyze_model_example_infinite_network_set(
    network_set: InfiniteHypergraphSet,
    calculated_properties: list[BaseNetworkProperty],
    save_directory: Path
) -> None:
    """Analyze the given infinite network set."""
    info('Infinite network set analysis started.')

    plt.rcParams["text.usetex"] = False

    summary = network_set.calc_network_summary(calculated_properties)

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


def create_infinite_network_plots(network: InfiniteHypergraph, save_path: Path) -> None:
    plot_network(network, False, save_path / 'infinite_network.png')
    plot_network(network, True, save_path / 'infinite_network_fixed_vertex_positions.png')
    plot_hypergraph(network, True, save_path / 'infinite_hypergraph_fixed_positions.png')
