#!/usr/bin/env python3
"""Helper functions for plotting."""

from enum import Enum, auto
from logging import warning
from typing import Any

import numpy as np
import numpy.typing as npt
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt

from distribution.approximation import DistributionApproximation
from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import create_fitting_parameters
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class PaddingSide(Enum):
    """Define which direction theoretical values should be padded."""

    NONE: int = auto()
    LEFT: int = auto()
    RIGHT: int = auto()
    BOTH: int = auto()


def approximate_and_plot_pdf(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    theoretical_distribution_type: TheoreticalDistribution.Type,
) -> None:
    """Plot the empirical distribution PDF and its approximation."""
    distribution_pair = DistributionApproximation(empirical_distribution, theoretical_distribution_type)
    default_fitting_parameters = create_fitting_parameters(theoretical_distribution_type)
    distribution_pair.fit(default_fitting_parameters)

    plot_distribution_pdf_pair(distribution_pair, axes)


def plot_distribution_pdf_pair(
    distribution_pair: DistributionApproximation,
    axes: plt.Axes
) -> None:
    """Plot the distribution and its approximation on a given axes."""
    # standardize the plots if the theoretical distribution is normal or stable
    if distribution_pair.type == TheoreticalDistribution.Type.POISSON:
        _plot_approximation(
            distribution_pair,
            EmpiricalDistribution.HistogramType.INTEGERS,
            PaddingSide.RIGHT,
            axes
        )
    elif distribution_pair.type == TheoreticalDistribution.Type.NORMAL:
        _plot_approximation_standardized(
            distribution_pair,
            EmpiricalDistribution.HistogramType.LINEAR,
            PaddingSide.BOTH,
            axes
        )
    elif distribution_pair.type == TheoreticalDistribution.Type.POWER_LAW:
        _plot_approximation_log(
            distribution_pair,
            EmpiricalDistribution.HistogramType.LOGARITHMIC,
            PaddingSide.NONE,
            axes
        )
    elif distribution_pair.type == TheoreticalDistribution.Type.STABLE:
        _plot_approximation_standardized(
            distribution_pair,
            EmpiricalDistribution.HistogramType.LINEAR,
            PaddingSide.BOTH,
            axes
        )


def _plot_approximation(
    distribution_pair: DistributionApproximation,
    histogram_type: EmpiricalDistribution.HistogramType,
    padding: PaddingSide,
    axes: plt.Axes
) -> None:

    if distribution_pair.empirical.valid:
        histogram, bins = distribution_pair.empirical.calc_histogram(histogram_type)
        _plot_histogram(histogram, bins, axes)
    else:
        warning('Empirical distribution to be plotted is invalid.')
        return

    if distribution_pair.theoretical.valid:
        domain_intersection = distribution_pair.theoretical.domain.intersect(distribution_pair.empirical.domain)
        theoretical_x_values = _get_theoretical_x_values(domain_intersection, padding)
        theoretical_pdf = distribution_pair.theoretical.pdf(theoretical_x_values)
        _plot_pdf(theoretical_x_values, theoretical_pdf, axes)
    else:
        theoretical_x_values, theoretical_pdf = np.empty((0,)), np.empty((0,))
        warning('Theoretical distribution to be plotted is invalid.')

    _set_linear_scale_limits(
        x_values=np.r_[bins, theoretical_x_values],  # merge all x coordinates to be plotted
        y_values=np.r_[histogram, theoretical_pdf],  # merge all y coordinates to be plotted
        axes=axes
    )

    test_result = distribution_pair.run_test()

    print_info([distribution_pair, test_result], axes)


def _plot_approximation_standardized(
    distribution_pair: DistributionApproximation,
    histogram_type: EmpiricalDistribution.HistogramType,
    padding: PaddingSide,
    axes: plt.Axes
) -> None:

    if distribution_pair.empirical.valid:
        mu = distribution_pair.empirical.mean
        std = distribution_pair.empirical.std_dev
        histogram, bins = distribution_pair.empirical.calc_histogram(histogram_type)
        standardized_bins, standardized_histogram = _standardize_coordinates(bins, histogram, mu, std)
        _plot_histogram(standardized_histogram, standardized_bins, axes)
    else:
        warning('Empirical distribution to be plotted is invalid.')
        return

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

    _set_linear_scale_limits(
        x_values=np.r_[standardized_bins, standardized_theoretical_x_values],  # merge all x coordinates to be plotted
        y_values=np.r_[standardized_histogram, standardized_theoretical_pdf],  # merge all y coordinates to be plotted
        axes=axes
    )

    test_result = distribution_pair.run_test()
    print_info([distribution_pair, test_result], axes)


def _plot_approximation_log(
    distribution_pair: DistributionApproximation,
    histogram_type: EmpiricalDistribution.HistogramType,
    padding: PaddingSide,
    axes: plt.Axes
) -> None:

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

    _set_logarithmic_scale_limits(
        x_values=np.r_[bins, theoretical_x_values],  # merge all x_values to be plotted
        y_values=np.r_[histogram, theoretical_pdf_to_plot],  # merge all y_values to be plotted
        axes=axes
    )

    test_result = distribution_pair.run_test()
    print_info([distribution_pair, test_result], axes)


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
) -> None:
    """Plot the empirical density with histogram estimation on the given axes with information."""
    histogram, bins = distribution.calc_histogram(histogram_type)
    _plot_histogram(histogram, bins, axes)
    _set_linear_scale_limits(x_values=bins, y_values=histogram, axes=axes)
    print_info([distribution], axes)


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


def plot_empirical_distribution_value_counts(distribution: EmpiricalDistribution, axes: plt.Axes) -> None:
    """Plot the empirical density with value counts on the given axes."""
    axes.set_xlabel('value')
    axes.set_ylabel('value counts')

    if not distribution.valid:
        warning('Distribution to be plotted was invalid.')
        return

    value_counts = distribution.calc_value_counts()
    plot_value_counts(value_counts, axes)


def plot_value_counts(value_counts: npt.NDArray[np.float_ | np.int_], axes: plt.Axes) -> None:
    """Plot a list of values on a given axes."""
    markerline, _, _ = axes.stem(value_counts[:, 0], value_counts[:, 1])
    plt.setp(markerline, markersize=5)
    axes.set_xlabel('Index')
    axes.set_ylabel('Value')


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


def print_info(objects_to_print: list[Any], axes: plt.Axes) -> None:
    """Print a text box to the given axes."""
    text = '\n'.join(map(str, objects_to_print))
    text_box = AnchoredText(text, frameon=True, loc='lower left', pad=0.5)
    plt.setp(text_box.patch, boxstyle='round', facecolor='wheat', alpha=0.5)
    axes.add_artist(text_box)


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
    axes.xaxis.get_label().get_text() + ' (log)'
    axes.xaxis.get_label().get_text() + ' (log)'
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
