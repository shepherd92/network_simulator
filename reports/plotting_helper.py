#!/usr/bin/env python3
"""Helper functions for plotting."""

from enum import Enum, auto
from logging import warning, debug
from operator import itemgetter
from pathlib import Path
from typing import Any

from gudhi.persistence_graphical_tools import (
    plot_persistence_barcode,
    plot_persistence_diagram,
)
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.collections import PolyCollection
import networkx as nx
import numpy as np
import numpy.typing as npt
from scipy.spatial import ConvexHull

from distribution.approximation import DistributionApproximation
from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from network.finite_network import FiniteNetwork
from tools.logging_helper import log_function_name


class PaddingSide(Enum):
    """Define which direction theoretical values should be padded."""

    NONE: int = auto()
    LEFT: int = auto()
    RIGHT: int = auto()
    BOTH: int = auto()


@log_function_name
def save_axes_as_separate_figure(file_path: Path, axes: plt.axes) -> None:
    """Save a particular axes as a separate figure."""
    figure = plt.figure()
    figure.axes.append(axes)
    figure.savefig(file_path)


@log_function_name
def plot_hypergraph(network: FiniteNetwork, determined_positions: bool, axes: plt.Axes) -> None:
    """Plot a hypergraph complex on the given axis."""
    if determined_positions:
        assert network.vertex_positions is not None
        assert network.interaction_positions is not None
        vertex_positions = network.vertex_positions
        interaction_positions = network.interaction_positions
    else:
        raise NotImplementedError('Determining positions is not implemented yet.')

    axes.scatter(
        [position[0] for position in vertex_positions.values()],
        [position[1] for position in vertex_positions.values()],
        c='blue', s=100
    )
    axes.scatter(
        [position[0] for position in interaction_positions.values()],
        [position[1] for position in interaction_positions.values()],
        c='red', s=100
    )

    # create the hypergraph network
    interactions = network.interactions
    num_of_interactions = len(interactions)
    hypergraph_edges = []
    for interaction_id, interaction in enumerate(network.interactions):
        hypergraph_edges.extend(list(zip(
            [interaction_id] * len(interaction),
            [vertex_id + num_of_interactions for vertex_id in interaction]
        )))

    # increment the vertex id-s by the number of interactions
    vertex_positions_shifted_id = {key + len(interactions): value for key, value in vertex_positions.items()}
    hypergraph_vertex_positions = interaction_positions | vertex_positions_shifted_id
    hypergraph = nx.Graph()
    hypergraph.add_nodes_from(hypergraph_vertex_positions.keys())
    hypergraph.add_edges_from(hypergraph_edges)

    nx.draw_networkx_edges(hypergraph, hypergraph_vertex_positions, ax=axes,
                           edge_color='black', width=1, alpha=1)
    axes.axhline(y=0., color='black', linestyle='-')
    axes.axhline(y=1., color='black', linestyle='-')
    axes.axvline(x=0., color='black', linestyle='-')
    axes.set_ylim(0., 1.)


@log_function_name
def plot_finite_network(network: FiniteNetwork, determined_vertex_positions: bool, axes: plt.Axes) -> None:
    """Plot a simplicial complex on the given axis."""
    color_map_name = 'plasma_r'  # viridis, plasma, inferno, magma, cividis
    if determined_vertex_positions:
        assert network.vertex_positions is not None
        vertex_positions = network.vertex_positions
    else:
        vertex_positions = _determine_node_positions(network.graph)
    interactions_to_plot = _determine_interactions_to_plot(network.interactions)
    debug(f'Number of interactions to plot: {len(interactions_to_plot)}')

    interaction_vertex_positions: list[list[tuple[float]]] = []
    for interaction in interactions_to_plot:
        this_interaction_vertex_positions = list(itemgetter(*interaction)(vertex_positions))
        if determined_vertex_positions:
            size_of_this_interaction = \
                max([coordinates[0] for coordinates in this_interaction_vertex_positions]) - \
                min([coordinates[0] for coordinates in this_interaction_vertex_positions])
            if size_of_this_interaction < 0.5:
                # exclude polygons which were created by wrap around torus effect
                interaction_vertex_positions.append(this_interaction_vertex_positions)
        else:
            interaction_vertex_positions.append(this_interaction_vertex_positions)

    convex_hulls = [
        ConvexHull(simplex_node_pos, qhull_options='QJ Pp')
        for simplex_node_pos in interaction_vertex_positions
    ]

    polygon_coordinates = [
        np.array([
            node_positions[node_index]
            for node_index in convex_hull.vertices
        ])
        for convex_hull, node_positions in zip(convex_hulls, interaction_vertex_positions)
    ]

    if len(interactions_to_plot) > 0:
        face_colors = _get_simplex_colors(interactions_to_plot, color_map_name)

        polygon_collection = PolyCollection(
            polygon_coordinates,
            facecolors=face_colors,
            edgecolors=('black',),
            linewidths=(0.0000001,)
        )

        axes.add_collection(polygon_collection)

    axes.set_xlim([
        min([coordinates[0] for coordinates in vertex_positions.values()]),
        max([coordinates[0] for coordinates in vertex_positions.values()]),
    ])
    axes.set_ylim([
        min([coordinates[1] for coordinates in vertex_positions.values()]),
        max([coordinates[1] for coordinates in vertex_positions.values()]),
    ])

    nx.draw_networkx_edges(network.graph, vertex_positions, ax=axes, edge_color='black', width=0.0001, alpha=1e-6)
    nx.draw_networkx_nodes(network.graph, vertex_positions, ax=axes, node_color='black', node_size=0.0001, alpha=1e-6)

    axes.set_axis_off()


@log_function_name
def plot_distribution_approximation(
    distribution_pair: DistributionApproximation,
    data_set_value: float,
    axes: plt.Axes,
) -> None:
    """Plot the distribution and its approximation on a given axes."""
    # standardize the plots if the theoretical distribution is normal or stable
    if distribution_pair.type == TheoreticalDistribution.Type.POISSON:
        plot_approximation_histogram(
            distribution_pair,
            EmpiricalDistribution.HistogramType.INTEGERS,
            data_set_value,
            PaddingSide.RIGHT,
            axes
        )
    elif distribution_pair.type == TheoreticalDistribution.Type.NORMAL:
        plot_approximation_histogram_standardized(
            distribution_pair,
            EmpiricalDistribution.HistogramType.LINEAR,
            data_set_value,
            PaddingSide.BOTH,
            axes
        )
    elif distribution_pair.type == TheoreticalDistribution.Type.POWER_LAW:
        plot_approximation_histogram_log(
            distribution_pair,
            EmpiricalDistribution.HistogramType.INTEGERS,
            data_set_value,
            PaddingSide.NONE,
            axes
        )
    elif distribution_pair.type == TheoreticalDistribution.Type.STABLE:
        plot_approximation_histogram_standardized(
            distribution_pair,
            EmpiricalDistribution.HistogramType.LINEAR,
            data_set_value,
            PaddingSide.BOTH,
            axes
        )


@log_function_name
def _determine_node_positions(graph: nx.Graph):

    if graph.number_of_nodes() <= 1000:
        graphviz_method = 'neato'  # best for a few nodes
    else:
        graphviz_method = 'sfdp'  # scalable method
    node_positions = nx.nx_agraph.graphviz_layout(graph, prog=graphviz_method)
    return node_positions


@log_function_name
def _determine_interactions_to_plot(interactions: list[list[int]]) -> list[list[int]]:
    # remove facets with dimension 0 and 1 (points and edges)
    interactions_to_plot = sorted([facet for facet in interactions if len(facet) - 1 > 1], key=len, reverse=False)
    return interactions_to_plot


@log_function_name
def _get_simplex_colors(simplices: list[list[int]], color_map_name: str):
    dimensions = np.array([len(simplex) - 1 for simplex in simplices])
    color_values_on_map = {dimension: (dimensions < dimension).mean() for dimension in np.unique(dimensions)}
    color_map = plt.get_cmap(color_map_name,)
    face_colors = color_map(np.vectorize(color_values_on_map.get)(dimensions))

    opacities = 1. / (dimensions - 1.)**2.
    face_colors[:, -1] = opacities
    return face_colors


@log_function_name
def plot_approximation_histogram(
    distribution_pair: DistributionApproximation,
    histogram_type: EmpiricalDistribution.HistogramType,
    data_set_value: float,
    padding: PaddingSide,
    axes: plt.Axes
) -> None:
    """Plot the histogram and the theoretical distribution on a given axes."""

    if distribution_pair.empirical.valid:
        histogram, bins = distribution_pair.empirical.calc_histogram(histogram_type)
        _plot_histogram(histogram, bins, axes)
    else:
        warning('Empirical distribution to be plotted is invalid.')
        return None

    if distribution_pair.theoretical.valid:
        domain_intersection = distribution_pair.theoretical.domain.intersect(distribution_pair.empirical.domain)
        theoretical_x_values = _get_theoretical_x_values(domain_intersection, padding)
        theoretical_pdf = distribution_pair.theoretical.pdf(theoretical_x_values)
        _plot_pdf(theoretical_x_values, theoretical_pdf, axes)
    else:
        theoretical_x_values, theoretical_pdf = np.empty((0,)), np.empty((0,))
        warning('Theoretical distribution to be plotted is invalid.')

    test_result = distribution_pair.run_test(data_set_value)
    _plot_point_value(data_set_value, axes)

    # merge all values to be plotted
    _set_linear_scale_limits(
        x_values=np.r_[bins, theoretical_x_values, data_set_value],
        y_values=np.r_[histogram, theoretical_pdf],
        axes=axes
    )

    print_info([distribution_pair, test_result], axes)


@log_function_name
def plot_approximation_histogram_standardized(
    distribution_pair: DistributionApproximation,
    histogram_type: EmpiricalDistribution.HistogramType,
    data_set_value: float,
    padding: PaddingSide,
    axes: plt.Axes
) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_]] | None:

    if distribution_pair.empirical.valid:
        mu = distribution_pair.empirical.mean
        std = distribution_pair.empirical.std_dev
        histogram, bins = distribution_pair.empirical.calc_histogram(histogram_type)
        standardized_bins, standardized_histogram = _standardize_coordinates(bins, histogram, mu, std)
        _plot_histogram(standardized_histogram, standardized_bins, axes)
    else:
        warning('Empirical distribution to be plotted is invalid.')
        return None

    if distribution_pair.theoretical.valid:
        domain_intersection = distribution_pair.theoretical.domain.intersect(distribution_pair.empirical.domain)
        theoretical_x_values = _get_theoretical_x_values(domain_intersection, padding)
        theoretical_pdf = distribution_pair.theoretical.pdf(theoretical_x_values)
        standardized_theoretical_x_values, standardized_theoretical_pdf = \
            _standardize_coordinates(theoretical_x_values, theoretical_pdf, mu, std)
        _plot_pdf(standardized_theoretical_x_values, standardized_theoretical_pdf, axes)
    else:
        standardized_theoretical_x_values, standardized_theoretical_pdf = np.empty((0,)), np.empty((0,))
        warning('Theoretical distribution to be plotted is invalid.')

    test_result = distribution_pair.run_test(data_set_value)
    _plot_point_value((data_set_value - mu) / std, axes)

    # merge all coordinates to be plotted
    _set_linear_scale_limits(
        x_values=np.r_[
            standardized_bins,
            standardized_theoretical_x_values,
            (data_set_value - mu) / std,
        ],
        y_values=np.r_[
            standardized_histogram,
            standardized_theoretical_pdf,
        ],
        axes=axes
    )

    print_info([distribution_pair, test_result], axes)

    values_plotted_empirical = np.c_[standardized_bins[:-1], standardized_histogram]
    values_plotted_theoretical = np.c_[standardized_theoretical_x_values, standardized_theoretical_pdf]

    return values_plotted_empirical, values_plotted_theoretical


@log_function_name
def plot_value_counts_log(
    value_counts: npt.NDArray[np.int_],
    axes: plt.Axes
) -> npt.NDArray[np.int_]:

    plot_value_counts(value_counts, axes)

    _set_logarithmic_scale_limits(
        x_values=value_counts[:, 0],
        y_values=value_counts[:, 1],
        axes=axes,
    )

    return value_counts


@log_function_name
def plot_approximation_value_counts_log(
    distribution_pair: DistributionApproximation,
    data_set_value: float,
    padding: PaddingSide,
    axes: plt.Axes
) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_]] | None:

    if distribution_pair.empirical.valid:
        value_counts = distribution_pair.empirical.calc_value_counts()
        plot_value_counts(value_counts, axes)
    else:
        value_counts = np.empty((0, 2))
        warning('Empirical distribution to be plotted is invalid.')

    if distribution_pair.theoretical.valid:
        domain_intersection = distribution_pair.theoretical.domain.intersect(distribution_pair.empirical.domain)

        theoretical_x_values = _get_theoretical_x_values(domain_intersection, padding)
        theoretical_pdf = distribution_pair.theoretical.pdf(theoretical_x_values)

        theoretical_integral = \
            distribution_pair.theoretical.cdf(np.array([domain_intersection.max_])) - \
            distribution_pair.theoretical.cdf(np.array([domain_intersection.min_]))

        empirical_integral = value_counts[
            (value_counts[:, 0] >= domain_intersection.min_) &
            (value_counts[:, 0] <= domain_intersection.max_)
        ][:, 1].sum()

        # correct for not identical domains: pdf-s should match in the domain intersection
        theoretical_pdf_to_plot = theoretical_pdf * empirical_integral / theoretical_integral
        _plot_pdf(theoretical_x_values, theoretical_pdf_to_plot, axes)
    else:
        theoretical_x_values, theoretical_pdf_to_plot = np.empty((0,)), np.empty((0,))
        warning('Theoretical distribution to be plotted is invalid.')

    test_result = distribution_pair.run_test(data_set_value)
    _plot_point_value(data_set_value, axes)

    # merge all values to be plotted
    _set_logarithmic_scale_limits(
        x_values=np.r_[
            value_counts[:, 0],
            theoretical_x_values,
            data_set_value,
        ],
        y_values=np.r_[
            value_counts[:, 1],
            theoretical_pdf_to_plot,
        ],
        axes=axes
    )

    print_info([distribution_pair, test_result], axes)

    values_plotted_empirical = value_counts
    values_plotted_theoretical = np.c_[theoretical_x_values, theoretical_pdf_to_plot]

    return values_plotted_empirical, values_plotted_theoretical


def plot_approximation_histogram_log(
    distribution_pair: DistributionApproximation,
    histogram_type: EmpiricalDistribution.HistogramType,
    data_set_value: float,
    padding: PaddingSide,
    axes: plt.Axes
) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_]] | None:

    if distribution_pair.empirical.valid:
        histogram, bins = distribution_pair.empirical.calc_histogram(histogram_type)
        _plot_histogram(histogram, bins, axes)
    else:
        histogram, bins = np.empty((0,)), np.empty((0,))
        warning('Empirical distribution to be plotted is invalid.')

    if distribution_pair.theoretical.valid:
        domain_intersection = distribution_pair.theoretical.domain.intersect(distribution_pair.empirical.domain)

        theoretical_x_values = _get_theoretical_x_values(domain_intersection, padding)
        theoretical_pdf = distribution_pair.theoretical.pdf(theoretical_x_values)

        theoretical_integral = \
            distribution_pair.theoretical.cdf(np.array([domain_intersection.max_])) - \
            distribution_pair.theoretical.cdf(np.array([domain_intersection.min_]))
        empirical_integral = \
            distribution_pair.empirical.cdf(np.array([domain_intersection.max_])) - \
            distribution_pair.empirical.cdf(np.array([domain_intersection.min_]))

        # correct for not identical domains: pdf-s should match in the domain intersection
        theoretical_pdf_to_plot = theoretical_pdf * empirical_integral / theoretical_integral
        _plot_pdf(theoretical_x_values, theoretical_pdf_to_plot, axes)
    else:
        theoretical_x_values, theoretical_pdf_to_plot = np.empty((0,)), np.empty((0,))
        warning('Theoretical distribution to be plotted is invalid.')

    test_result = distribution_pair.run_test(data_set_value)
    _plot_point_value(data_set_value, axes)

    # merge all values to be plotted
    _set_logarithmic_scale_limits(
        x_values=np.r_[bins, theoretical_x_values, data_set_value],
        y_values=np.r_[histogram, theoretical_pdf_to_plot],
        axes=axes
    )

    print_info([distribution_pair, test_result], axes)

    values_plotted_empirical = np.c_[bins[:-1], histogram]
    values_plotted_theoretical = np.c_[theoretical_x_values, theoretical_pdf_to_plot]

    return values_plotted_empirical, values_plotted_theoretical


def _standardize_coordinates(
    x_values: npt.NDArray[np.float_],
    pdf: npt.NDArray[np.float_],
    location: np.float_,
    scale: np.float_
) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_]]:

    standardized_x_values = (x_values - location) / scale
    standardized_pdf = pdf * scale  # if the domain is rescaled, the values of the pdf also change

    return standardized_x_values, standardized_pdf


def plot_empirical_distribution_histogram_with_info(
    distribution: EmpiricalDistribution,
    histogram_type: EmpiricalDistribution.HistogramType,
    axes: plt.Axes,
) -> npt.NDArray[np.float_]:
    """Plot the empirical density with histogram estimation on the given axes with information."""
    histogram, bins = distribution.calc_histogram(histogram_type)
    _plot_histogram(histogram, bins, axes)
    _set_linear_scale_limits(x_values=bins, y_values=histogram, axes=axes)
    print_info([distribution], axes)

    values_plotted = np.c_[bins[:-1], histogram]
    return values_plotted


def _plot_pdf(x_values: npt.NDArray[np.float_], pdf: npt.NDArray[np.float_], axes: plt.Axes) -> None:
    axes.plot(x_values, pdf, color='red')


def _plot_histogram(
    histogram: npt.NDArray[np.float_],
    bins: npt.NDArray[np.float_],
    axes: plt.Axes
) -> None:
    """Plot the empirical density with histogram estimation on the given axes."""
    axes.set_xlabel('value')
    axes.set_ylabel('density')

    centers = (bins[:-1] + bins[1:]) / 2
    widths = 0.9 * np.diff(bins)

    axes.bar(centers, histogram, align='center', width=widths)


def _plot_point_value(value: float, axes: plt.Axes) -> None:
    """Plot a vertical line to the point value on the given axes."""
    axes.axvline(x=value, color='green')


def plot_empirical_distribution_value_counts(
    distribution: EmpiricalDistribution,
    axes: plt.Axes
) -> None:
    """Plot the empirical density with value counts on the given axes."""
    axes.set_xlabel('value')
    axes.set_ylabel('value counts')

    if not distribution.valid:
        warning('Distribution to be plotted was invalid.')
        return

    value_counts = distribution.calc_value_counts()
    plot_value_counts(value_counts, axes)


def plot_value_counts(
    value_counts: npt.NDArray[np.float_ | np.int_],
    axes: plt.Axes
) -> npt.NDArray[np.float_]:
    """Plot a list of values on a given axes."""
    # markerline, _, _ = axes.stem(value_counts[:, 0], value_counts[:, 1])
    # plt.setp(markerline, markersize=5)
    axes.scatter(value_counts[:, 0], value_counts[:, 1])
    _set_linear_scale_limits(
        x_values=value_counts[:, 0],
        y_values=value_counts[:, 1],
        axes=axes
    )
    axes.set_xlabel('Index')
    axes.set_ylabel('Value')
    return value_counts


def plot_probability_plot(distribution_pair: DistributionApproximation, axes: plt.Axes) -> None:
    """Plot the probability plot on a given axes."""
    probability_plot_points = distribution_pair.generate_probability_plot_points()
    axes.scatter(probability_plot_points[:, 0], probability_plot_points[:, 1])
    axes.plot([0, 1], [0, 1], color='red')

    test_results = distribution_pair.run_test()
    print_info([test_results], axes)

    axes.set_xlim([0, 1])
    axes.set_ylim([0, 1])
    axes.set_xlabel('theoretical')
    axes.set_ylabel('empirical')


def plot_qq_plot(distribution_pair: DistributionApproximation, axes: plt.Axes) -> None:
    """Plot the Q-Q plot on a given axes."""
    qq_plot_points = distribution_pair.generate_qq_plot_points()
    axes.scatter(qq_plot_points[:, 0], qq_plot_points[:, 1])

    test_results = distribution_pair.run_test()
    print_info([test_results], axes)

    axes.set_xlabel('theoretical quantiles')
    axes.set_ylabel('empirical quantiles')
    limits = [
        np.min([axes.get_xlim(), axes.get_ylim()]),  # min of both axes
        np.max([axes.get_xlim(), axes.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    axes.plot(limits, limits, color='red')

    axes.set_aspect('equal')
    axes.set_xlim(limits)
    axes.set_ylim(limits)


def plot_persistence_diagram_(persistence, axes: plt.Axes):
    """Plot the persistence diagram and return the plotted values."""
    plot_persistence_diagram(persistence, axes=axes)
    return persistence


def plot_persistence_barcode_(persistence, axes: plt.Axes):
    """Plot the persistence barcode and return the plotted values."""
    plot_persistence_barcode(persistence, axes=axes)
    return persistence


def print_info(objects_to_print: list[Any], axes: plt.Axes) -> None:
    """Print a text box to the given axes."""
    text = '\n'.join(map(str, objects_to_print))
    text_box = AnchoredText(text, frameon=True, loc='lower left', pad=0.5)
    plt.setp(text_box.patch, boxstyle='round', facecolor='wheat', alpha=0.5)
    axes.add_artist(text_box)


def check_calculated(plotter_function):
    """Check if the data to be plotted is calculated.

    Use as a decorator.
    """
    def wrapper(data_to_plot: Any, axes: plt.Axes, save_directory: Path):
        if data_to_plot is None:
            print_not_calculated(axes)
        else:
            plotter_function(data_to_plot, axes, save_directory)
    return wrapper


def print_not_calculated(axes: plt.Axes) -> None:
    """Print the text not calculated to the given axes."""
    properties = dict(fontweight='bold', fontsize=20)
    text_box = AnchoredText('NOT CALCULATED', frameon=True, loc='center', pad=0.5, prop=properties)
    plt.setp(text_box.patch, boxstyle='round', facecolor='wheat', alpha=0.5)
    axes.add_artist(text_box)


def _get_theoretical_x_values(
    plot_domain: Distribution.Domain,
    padding_sides: PaddingSide
) -> npt.NDArray[np.float_]:

    NUM_OF_POINTS = 100
    PADDING = 0.2

    if padding_sides == PaddingSide.NONE:
        # add padding to the right
        x_values = np.linspace(
            plot_domain.min_,
            plot_domain.max_,
            NUM_OF_POINTS,
            endpoint=True
        )
    elif padding_sides == PaddingSide.LEFT:
        # add padding to the right
        x_values = np.linspace(
            plot_domain.min_ - PADDING * plot_domain.length,
            plot_domain.max_,
            NUM_OF_POINTS,
            endpoint=True
        )
    elif padding_sides == PaddingSide.RIGHT:
        # add padding to the right
        x_values = np.linspace(
            plot_domain.min_,
            plot_domain.max_ + PADDING * plot_domain.length,
            NUM_OF_POINTS,
            endpoint=True
        )
    elif padding_sides == PaddingSide.BOTH:
        # add padding to the right
        x_values = np.linspace(
            plot_domain.min_ - PADDING * plot_domain.length,
            plot_domain.max_ + PADDING * plot_domain.length,
            NUM_OF_POINTS,
            endpoint=True
        )

    return x_values


def _set_linear_scale_limits(
    x_values: npt.NDArray[np.float_ | np.int_],
    y_values: npt.NDArray[np.float_ | np.int_],
    axes: plt.Axes
) -> None:
    x_min, x_max = _calc_linear_scale_plot_limits(x_values)
    y_min, y_max = _calc_linear_scale_plot_limits(y_values)

    axes.set_xlim([x_min, x_max])
    axes.set_ylim([y_min, y_max])


def _set_logarithmic_scale_limits(
    x_values: npt.NDArray[np.float_ | np.int_],
    y_values: npt.NDArray[np.float_ | np.int_],
    axes: plt.Axes
) -> None:
    axes.set_xlabel(axes.xaxis.get_label().get_text() + ' (log)')
    axes.set_ylabel(axes.xaxis.get_label().get_text() + ' (log)')
    axes.set_xscale('log')
    axes.set_yscale('log')

    x_min, x_max = _calc_log_scale_plot_limits(x_values)
    y_min, y_max = _calc_log_scale_plot_limits(y_values)
    axes.set_xlim([x_min, x_max])
    axes.set_ylim([y_min, y_max])


def _calc_log_scale_plot_limits(values_single_axis: npt.NDArray[np.float_ | np.int_]) -> tuple[float, float]:

    padding = 0.1

    finite_values = values_single_axis[np.isfinite(values_single_axis)]
    positive_finite_values: npt.NDArray[np.float_] = finite_values[finite_values > 0.]

    if len(positive_finite_values) == 0:
        # no finite values
        lower_limit = 1.
        upper_limit = 10.
    elif len(np.unique(positive_finite_values)) == 1:
        # only a single finite value
        lower_limit = positive_finite_values[0] * (1 - padding)
        upper_limit = positive_finite_values[0] * (1 + padding)
    else:
        values_min: float = positive_finite_values.min()
        values_max: float = positive_finite_values.max()
        lower_limit = values_min**(1 + padding) / values_max**padding
        upper_limit = values_max**(1 + padding) / values_min**padding

    assert not np.isnan(lower_limit) and not np.isnan(upper_limit)

    return lower_limit, upper_limit


def _calc_linear_scale_plot_limits(values_single_axis: npt.NDArray[np.float_ | np.int_]) -> tuple[float, float]:

    padding = 0.1

    finite_values: npt.NDArray[np.float_] = values_single_axis[np.isfinite(values_single_axis)]
    if len(finite_values) == 0:
        # no finite values
        lower_limit = 0.
        upper_limit = 1.
    elif len(np.unique(finite_values)) == 1:
        # only a single finite value
        lower_limit = finite_values[0] - 0.1
        upper_limit = finite_values[0] + 0.1
    else:
        values_min: float = finite_values.min()
        values_max: float = finite_values.max()
        lower_limit = values_min - padding * (values_max - values_min)
        upper_limit = values_max + padding * (values_max - values_min)

    assert not np.isnan(lower_limit) and not np.isnan(upper_limit)

    return lower_limit, upper_limit


def plot_giant_component(network: FiniteNetwork, save_path: Path):
    """Plot the giant component of the network."""
    plt.rcParams["text.usetex"] = False

    debug('Plotting simplicial complex of giant component started.')
    figure, axes = plt.subplots(1, 1)
    network_to_plot = network.get_component(0)
    plot_finite_network(network_to_plot, False, axes)
    figure.savefig(save_path)
    figure.clf()
    debug('Plotting simplicial complex of giant component finished.')


def plot_network_determined_positions(network: FiniteNetwork, save_path: Path):
    """Plot the entire network with predetermined vertex positions."""
    plt.rcParams["text.usetex"] = False

    debug('Plotting simplicial complex with fixed positions started.')
    figure, axes = plt.subplots(1, 1, figsize=(50, 50))
    plot_finite_network(network, True, axes)
    figure.savefig(save_path)
    figure.clf()
    debug('Plotting simplicial complex with fixed positions finished.')


def plot_hypergraph_determined_positions(network: FiniteNetwork, save_path: Path):
    """Plot the entire network with predetermined vertex positions."""
    plt.rcParams["text.usetex"] = False

    debug('Plotting hypergraph fixed positions started.')
    figure, axes = plt.subplots(1, 1, figsize=(50, 50))
    plot_hypergraph(network, True, axes)
    figure.savefig(save_path)
    figure.clf()
    debug('Plotting hypergraph fixed positions finished.')
