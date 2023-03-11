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
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution


class PlotScale(Enum):
    """Define the scaling of a plot axis."""

    LINEAR: int = auto()
    LOGARITHMIC: int = auto()


def approximate_and_plot_pdf(
    empirical_distribution: EmpiricalDistribution,
    axes: plt.Axes,
    theoretical_distribution_type: TheoreticalDistribution.Type,
) -> None:
    """Plot the empirical distribution PDF and its approximation."""
    distribution_pair = DistributionApproximation(empirical_distribution, theoretical_distribution_type)
    distribution_pair.fit()
    plot_distribution_pdf_pair(distribution_pair, axes)


def plot_distribution_pdf_pair(distribution_pair: DistributionApproximation, axes: plt.Axes) -> None:
    """Plot the distribution and its approximation on a given axes."""
    plot_empirical_distribution_pdf(distribution_pair.empirical, axes)

    x_values = distribution_pair.empirical.natural_x_values
    empirical_pdf = distribution_pair.empirical.pdf(x_values)

    if len(x_values) == 0 or len(empirical_pdf[np.isfinite(empirical_pdf)]) == 0:
        return

    domain_intersection = distribution_pair.theoretical.domain.intersect(distribution_pair.empirical.domain)
    theoretical_x_values: npt.NDArray[np.float_] = np.linspace(
        domain_intersection.min_,
        domain_intersection.max_,
        100,
        endpoint=True
    )

    theoretical_pdf = distribution_pair.theoretical.pdf(theoretical_x_values)

    theoretical_integral = \
        distribution_pair.theoretical.cdf(np.array([domain_intersection.max_])) - \
        distribution_pair.theoretical.cdf(np.array([domain_intersection.min_]))
    empirical_integral = \
        distribution_pair.empirical.cdf(np.array([domain_intersection.max_])) - \
        distribution_pair.empirical.cdf(np.array([domain_intersection.min_]))

    # empirical_integral = np.trapz(empirical_pdf, x_values)
    # theoretical_integral = np.trapz(theoretical_pdf, theoretical_x_values)

    theoretical_pdf_to_plot = theoretical_pdf * empirical_integral / theoretical_integral

    axes.plot(theoretical_x_values, theoretical_pdf_to_plot, color='red')

    if isinstance(distribution_pair.theoretical, PowerLawDistribution):
        axes.set_xlabel('value (log)')
        axes.set_ylabel('density (log)')
        axes.set_xscale('log')
        axes.set_yscale('log')
        x_min, x_max = _calc_log_scale_plot_limits(x_values)
        y_min, y_max = _calc_log_scale_plot_limits(empirical_pdf)
    else:
        axes.set_xlabel('value')
        axes.set_ylabel('density')
        x_min, x_max = _calc_linear_scale_plot_limits(x_values)
        y_min, y_max = _calc_linear_scale_plot_limits(empirical_pdf)

    axes.set_xlim([x_min, x_max])
    axes.set_ylim([y_min, y_max])

    test_result = distribution_pair.run_test()

    print_info([distribution_pair, test_result], axes)


def plot_empirical_distribution_pdf_with_info(distribution: EmpiricalDistribution, axes: plt.Axes) -> None:
    """Plot the distribution and its approximation on a given axes."""
    plot_empirical_distribution_pdf(distribution, axes)
    print_info([distribution], axes)


def plot_empirical_distribution_pdf(distribution: EmpiricalDistribution, axes: plt.Axes) -> None:
    """Plot the distribution and its approximation on a given axes."""
    axes.set_xlabel('value')
    axes.set_ylabel('density')

    if not distribution.valid:
        warning('Distribution to be plotted was invalid.')
        return

    x_values = distribution.natural_x_values
    pdf = distribution.pdf(x_values)

    if len(x_values) == 0 or len(pdf[np.isfinite(pdf)]) == 0:
        warning('Distribution to be plotted has no finite values.')
        return

    axes.plot(x_values, pdf)

    x_min, x_max = _calc_linear_scale_plot_limits(x_values)
    y_min, y_max = _calc_linear_scale_plot_limits(pdf)
    axes.set_xlim([x_min, x_max])
    axes.set_ylim([y_min, y_max])


def plot_value_list(value_list: list[int | float], axes: plt.Axes) -> None:
    """Plot a list of values on a given axes."""
    markerline, _, _ = axes.stem(range(len(value_list)), value_list)
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

    axes.set_xlabel('theoretical')
    axes.set_ylabel('empirical')
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
