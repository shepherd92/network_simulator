#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from logging import info
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import (
    create_fitting_parameters_power_law_data_set,
)
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from network.finite_network import FiniteNetwork
from network.property import BaseNetworkProperty
from reports.plotting_helper import (
    PaddingSide,
    check_calculated,
    plot_approximation_value_counts_log,
    plot_giant_component,
    plot_persistence_barcode_,
    plot_persistence_diagram_,
    plot_value_counts,
    plot_value_counts_log,
    print_not_calculated,
)
from tools.logging_helper import log_function_name


@log_function_name
def analyze_data_set_network(
    network: FiniteNetwork,
    calculated_properties: list[BaseNetworkProperty.Type],
    plot: bool,
    save_directory: Path,
) -> None:
    """Generate exploratory analysis of the provided finite network."""
    info('Finite network analysis started.')

    plt.rcParams["text.usetex"] = False

    if plot:
        plot_giant_component(network, save_directory / 'network.png')

    summary = network.calc_network_summary(calculated_properties)

    axes_grid_height = 4
    axes_grid_width = 3

    network.save_info(save_directory / 'network_info.csv')

    figure = plt.figure('Network Analysis', figsize=(axes_grid_width * 10, axes_grid_height * 10))
    axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)
    subfigure_row_index = 0

    if BaseNetworkProperty.Type.EDGES in calculated_properties:
        edges = summary.get(BaseNetworkProperty.Type.EDGES)
        np.savetxt(save_directory / 'edges.csv', edges, delimiter=',')

    _report_degree_distribution(
        summary.get(BaseNetworkProperty.Type.DEGREE_DISTRIBUTION),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory
    )
    _report_in_degree_distribution(
        summary.get(BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory
    )
    _report_out_degree_distribution(
        summary.get(BaseNetworkProperty.Type.OUT_DEGREE_DISTRIBUTION),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )

    subfigure_row_index += 1

    _report_vertices_by_component(
        summary.get(BaseNetworkProperty.Type.VERTICES_BY_COMPONENT),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory
    )
    _report_ho_1_degree_distribution(
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory
    )
    _report_ho_2_degree_distribution(
        summary.get(BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )

    subfigure_row_index += 1

    _report_simplex_dimension_distribution(
        summary.get(BaseNetworkProperty.Type.SIMPLEX_DIMENSION_DISTRIBUTION),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory
    )
    # _report_facet_dimension_distribution(
    #     summary.get(BaseNetworkProperty.Type.FACET_DIMENSION_DISTRIBUTION),
    #     figure.add_subplot(axes_grid[subfigure_row_index, 1]),
    #     save_directory
    # )
    _report_interaction_degree_distribution(
        summary.get(BaseNetworkProperty.Type.VERTEX_INTERACTION_DEGREE_DISTRIBUTION),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory
    )
    _report_interaction_dimension_distribution(
        summary.get(BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )

    subfigure_row_index += 1

    _report_betti_numbers(
        summary.get(BaseNetworkProperty.Type.BETTI_NUMBERS),
        figure.add_subplot(axes_grid[subfigure_row_index, 0]),
        save_directory
    )
    _report_persistence_diagram(
        summary.get(BaseNetworkProperty.Type.PERSISTENCE),
        figure.add_subplot(axes_grid[subfigure_row_index, 1]),
        save_directory
    )
    _report_betti_number_1_by_component(
        summary.get(BaseNetworkProperty.Type.BETTI_NUMBERS_BY_COMPONENT),
        figure.add_subplot(axes_grid[subfigure_row_index, 2]),
        save_directory
    )

    subfigure_row_index += 1

    figure.tight_layout()
    figure.savefig(save_directory / 'whole_report.png')
    figure.clf()

    info('Finite network analysis finished.')


@check_calculated
def _report_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'total_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)
    axes.set_title('Total Degree Distribution')


@check_calculated
def _report_in_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'in_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)
    axes.set_title('In Degree Distribution')


@check_calculated
def _report_out_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'out_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.RIGHT, axes)
    axes.set_title('Out Degree Distribution')


@check_calculated
def _report_ho_1_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'ho_degree_distribution_1')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)
    axes.set_title('Higher-Order Degree Distribution - Dimension 1')


@check_calculated
def _report_ho_2_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'ho_degree_distribution_2')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)
    axes.set_title('Higher-Order Degree Distribution - Dimension 2')


@check_calculated
def _report_ho_3_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'ho_degree_distribution_3')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)
    axes.set_title('Higher-Order Degree Distribution - Dimension 3')


@check_calculated
def _report_simplex_dimension_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'simplex_dimension_distribution')
    plot_value_counts(empirical_distribution.calc_value_counts(), axes)
    axes.set_title('Simplex dimension distribution')


@check_calculated
def _report_facet_dimension_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'facet_dimension_distribution')
    plot_value_counts(empirical_distribution.calc_value_counts(), axes)
    axes.set_title('Facet dimension distribution')


@check_calculated
def _report_interaction_dimension_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'interaction_dimension_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)
    axes.set_title('Interaction dimension distribution')


@check_calculated
def _report_interaction_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'interaction_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)
    axes.set_title('Interaction degree distribution')


@check_calculated
def _report_betti_numbers(
    betti_numbers: list[int],
    axes: plt.Axes,
    save_directory: Path
) -> None:
    axes.set_title('Betti Numbers')
    betti_numbers = np.c_[np.arange(len(betti_numbers)), np.array(betti_numbers)]
    plot_value_counts(betti_numbers, axes)
    pd.DataFrame(
        betti_numbers,
        columns=['dimension', 'betti_number'],
        dtype=np.int32,
    ).to_csv(save_directory / 'betti_numbers.csv', index=False)


@check_calculated
def _report_betti_number_1_by_component(
    betti_numbers_by_component: list[list[int]],
    axes: plt.Axes,
    save_directory: Path
) -> None:
    axes.set_title('Betti Number 1 by Component')

    betti_numbers_1 = np.array([
        betti_numbers_in_component[1]
        for betti_numbers_in_component in betti_numbers_by_component
    ])

    values_to_plot = np.c_[np.arange(len(betti_numbers_1)), betti_numbers_1,]
    plot_value_counts_log(values_to_plot, axes)
    data_frame = pd.DataFrame(
        betti_numbers_by_component,
        columns=list(range(betti_numbers_by_component.shape[1])),
        dtype=np.int32,
    )
    data_frame.index.name = 'component_index'
    data_frame.to_csv(save_directory / 'betti_numbers_by_component.csv')


@check_calculated
def _report_vertices_by_component(
    vertices_by_component: npt.NDArray[np.int_],
    axes: plt.Axes,
    save_directory: Path
) -> None:
    axes.set_title('Vertices by Component')

    values_to_plot = np.c_[
        np.arange(1, len(vertices_by_component) + 1),
        vertices_by_component,
    ]
    plot_value_counts_log(values_to_plot, axes)
    pd.DataFrame(
        values_to_plot,
        columns=['component_index', 'num_of_vertices'],
        dtype=np.int32,
    ).to_csv(save_directory / 'vertices_by_component.csv', index=False)


@check_calculated
def _report_persistence_diagram(
    persistence: tuple[tuple[int, tuple[int, int]], ...],
    axes: plt.Axes,
    save_directory: Path
) -> None:
    axes.set_title('Persistence Diagram')
    plot_persistence_diagram_(persistence, axes)
    _save_persistence(persistence, save_directory)


@check_calculated
def _report_persistence_barcode(
    persistence: tuple[tuple[int, tuple[int, int]], ...],
    axes: plt.Axes,
    save_directory: Path
) -> None:
    axes.set_title('Persistence Barcode')
    plot_persistence_barcode_(persistence, axes)
    _save_persistence(persistence, save_directory)


def _save_persistence(persistence: tuple[tuple[int, tuple[int, int]], ...], save_directory: Path) -> None:
    pd.DataFrame(
        [[dimension, birth, death] for dimension, (birth, death) in persistence],
        columns=['dimension', 'birth', 'death']
    ).to_csv(save_directory / 'persistence.csv', index=False)


def _plot_base_property(
    property_value: Any,
    title: str,
    plotter: Callable[[Any, plt.Axes], None],
    axes: plt.Axes
) -> None:

    axes.set_title(title)

    if property_value is None:
        print_not_calculated(axes)
        return None

    plotter(property_value, axes=axes)
