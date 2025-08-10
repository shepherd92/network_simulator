#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from logging import info
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from config_files.distribution_fitting_params import POWER_LAW_FITTING_MINIMUM_VALUE_DATA
from network.finite_network import FiniteNetwork
from network.property import BaseNetworkProperty
from reports.finite_network_analysis import (
    report_vertex_edge_degree_distribution,
    report_in_degree_distribution,
    report_out_degree_distribution,
    report_vertices_by_component,
    report_edge_triangle_degree_distribution,
    report_facet_dimension_distribution,
    report_triangle_tetrahedron_degree_distribution,
    report_simplex_dimension_distribution,
    report_vertex_interaction_degree_distribution,
    report_edge_interaction_degree_distribution,
    report_interaction_dimension_distribution,
    report_betti_numbers,
    report_persistence_diagram,
    # report_betti_number_1_by_component,
)
from reports.plotting_helper import plot_network
from tools.logging_helper import log_function_name


@log_function_name
def analyze_data_set_network(
    network: FiniteNetwork,
    calculated_properties: list[BaseNetworkProperty],
    plot: bool,
    save_directory: Path,
) -> None:
    """Generate exploratory analysis of the provided finite network."""
    info('Finite network analysis started.')

    plt.rcParams["text.usetex"] = False

    if plot:
        plot_network(network.get_component(0), False, save_directory / 'network.png')

    summary = network.calc_network_summary(calculated_properties)

    axes_grid_height = 5
    axes_grid_width = 3

    network.save_info(save_directory / 'network_info.csv')

    figure = plt.figure('Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
    axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
    subfigure_row_index = 0

    report_in_degree_distribution(
        summary.get(BaseNetworkProperty.in_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )
    report_out_degree_distribution(
        summary.get(BaseNetworkProperty.out_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )

    subfigure_row_index += 1

    report_vertex_edge_degree_distribution(
        summary.get(BaseNetworkProperty.vertex_edge_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )
    report_edge_triangle_degree_distribution(
        summary.get(BaseNetworkProperty.edge_triangle_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )
    report_triangle_tetrahedron_degree_distribution(
        summary.get(BaseNetworkProperty.triangle_tetrahedra_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )

    subfigure_row_index += 1

    report_vertex_interaction_degree_distribution(
        summary.get(BaseNetworkProperty.vertex_interaction_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )
    report_edge_interaction_degree_distribution(
        summary.get(BaseNetworkProperty.edge_interaction_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )
    report_interaction_dimension_distribution(
        summary.get(BaseNetworkProperty.interaction_vertex_degree_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )

    subfigure_row_index += 1

    report_simplex_dimension_distribution(
        summary.get(BaseNetworkProperty.simplex_dimension_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    )
    report_facet_dimension_distribution(
        summary.get(BaseNetworkProperty.facet_dimension_distribution),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory,
        power_law_fitting_minimum_value=POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
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
    figure.savefig(save_directory / 'whole_report.png')
    figure.clf()

    info('Finite network analysis finished.')
