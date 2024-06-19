from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import (
    create_fitting_parameters_poisson,
    create_power_law_fitting_parameters,
)
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from reports.plotting_helper import (
    PaddingSide,
    check_calculated,
    plot_approximation_value_counts_log,
    plot_persistence_barcode_,
    plot_persistence_diagram_,
    plot_value_counts,
    plot_value_counts_log,
)


@check_calculated('Vertex--edge degree distribution')
def report_vertex_edge_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    **kwargs,
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(kwargs['power_law_fitting_minimum_value'])
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'total_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


@check_calculated('In degree distribution')
def report_in_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    **kwargs,
) -> None:
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW,
    )
    fitting_parameters = create_power_law_fitting_parameters(kwargs['power_law_fitting_minimum_value'])
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'in_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


@check_calculated('Out degree distribution')
def report_out_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    **_,
) -> None:
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POISSON,
    )
    fitting_parameters = create_fitting_parameters_poisson()
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'out_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.RIGHT, axes)


@check_calculated('Edge--triangle degree distribution')
def report_edge_triangle_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    power_law_fitting_minimum_value: float,
) -> None:

    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(power_law_fitting_minimum_value)
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'ho_degree_distribution_1')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


@check_calculated('Triangle--tetrahedron degree distribution')
def report_triangle_tetrahedron_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    power_law_fitting_minimum_value: float,
) -> None:
    """Report triangle--tetrahedron degree distribution."""
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(power_law_fitting_minimum_value)
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'ho_degree_distribution_2')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


@check_calculated('Simplex dimension distribution')
def report_simplex_dimension_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    power_law_fitting_minimum_value: float,
) -> None:
    """Report simplex dimension distribution."""
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(power_law_fitting_minimum_value)
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'simplex_dimension_distribution')
    plot_value_counts(empirical_distribution.calc_value_counts(), axes)


@check_calculated('Facet dimension distribution')
def report_facet_dimension_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    power_law_fitting_minimum_value: float,
) -> None:
    """Report facet dimension distribution."""
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(power_law_fitting_minimum_value)
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'facet_dimension_distribution')
    plot_value_counts(empirical_distribution.calc_value_counts(), axes)


@check_calculated('Vertex--interaction degree distribution')
def report_vertex_interaction_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    power_law_fitting_minimum_value: float,
) -> None:
    """Report vertex--interaction degree distribution."""
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(power_law_fitting_minimum_value)
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'vertex_interaction_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


@check_calculated('Edge--interaction degree distribution')
def report_edge_interaction_degree_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    power_law_fitting_minimum_value: float,
) -> None:
    """Report edge--interaction degree distribution."""
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(power_law_fitting_minimum_value)
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'edge_interaction_degree_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


@check_calculated('Interaction--vertex degree distribution')
def report_interaction_dimension_distribution(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    save_directory: Path,
    power_law_fitting_minimum_value: float,
) -> None:
    """Report interaction--vertex degree distribution."""
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_power_law_fitting_parameters(power_law_fitting_minimum_value)
    approximation.fit(fitting_parameters)

    approximation.save(save_directory / 'interaction_dimension_distribution')
    plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


@check_calculated('Betti numbers')
def report_betti_numbers(
    betti_numbers: list[int],
    axes: plt.Axes,
    save_directory: Path,
    **_,
) -> None:
    """Report Betti numbers."""
    betti_numbers_array = np.c_[np.arange(len(betti_numbers)), np.array(betti_numbers)]
    plot_value_counts(betti_numbers_array, axes)
    pd.DataFrame(
        betti_numbers_array,
        columns=['dimension', 'betti_number'],
        dtype=np.int32,
    ).to_csv(save_directory / 'betti_numbers.csv', index=False)


@check_calculated('Betti number 1 by component')
def report_betti_number_1_by_component(
    betti_numbers_by_component: list[list[int]],
    axes: plt.Axes,
    save_directory: Path,
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
    data_frame.to_csv(save_directory / 'betti_numbers_by_component.csv')


@check_calculated('Vertices by component')
def report_vertices_by_component(
    vertices_by_component: npt.NDArray[np.int_],
    axes: plt.Axes,
    save_directory: Path,
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
    ).to_csv(save_directory / 'vertices_by_component.csv', index=False)


@check_calculated('Persistence diagram')
def report_persistence_diagram(
    persistence: list[list[tuple[float, float]]],
    axes: plt.Axes,
    save_directory: Path,
    **_,
) -> None:
    """Report persistence diagram."""
    plot_persistence_diagram_(persistence, axes)
    _save_persistence(persistence, save_directory)


@check_calculated('Persistence barcode')
def report_persistence_barcode(
    persistence: list[list[tuple[float, float]]],
    axes: plt.Axes,
    save_directory: Path,
    **_,
) -> None:
    """Report persistence barcode."""
    plot_persistence_barcode_(persistence, axes)
    _save_persistence(persistence, save_directory)


def _save_persistence(persistence: list[list[tuple[float, float]]], save_directory: Path) -> None:
    persistence_df = pd.DataFrame([
            [dimension, birth, death]
            for dimension, intervals in enumerate(persistence)
            for birth, death in intervals
        ],
        columns=['dimension', 'birth', 'death']
    ).astype(int)

    persistence_df = persistence_df.groupby(persistence_df.columns.tolist(), as_index=False).size()

    persistence_df.to_csv(save_directory / 'persistence.csv', index=False)
