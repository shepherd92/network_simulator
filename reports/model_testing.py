#!/usr/bin/env python3
"""This module is responsible for generating a report for a model."""

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from network.property import ScalarNetworkPropertyReport
from reports.plotting_helper import plot_distribution_approximation, plot_probability_plot, plot_qq_plot


def create_model_test_report(scalar_property_reports: list[ScalarNetworkPropertyReport]) -> Figure:
    """Create a figure with plots of the given scalar property report."""
    axes_grid_height = len(scalar_property_reports)
    axes_grid_width = 2

    figure = plt.figure('Model test report', figsize=(axes_grid_width * 10, axes_grid_height * 10))
    axes_grid = figure.add_gridspec(axes_grid_height, axes_grid_width)

    for index, property_report in enumerate(scalar_property_reports):
        _plot_property_report(property_report, figure.add_subplot(axes_grid[index, 0]))
        _plot_qq_plot(property_report, figure.add_subplot(axes_grid[index, 1]))
        # _plot_probability_plot(property_report, figure.add_subplot(axes_grid[index, 2]))

    figure.tight_layout()

    return figure


def _plot_property_report(property_report: ScalarNetworkPropertyReport, axes: plt.Axes) -> None:
    """Plot the distribution and its approximation on a given axes."""
    axes.set_title(f'{property_report.params.name} distribution')

    distribution_pair = property_report.distributions
    plot_distribution_approximation(distribution_pair, property_report.data_point, axes)


def _plot_probability_plot(property_report: ScalarNetworkPropertyReport, axes: plt.Axes) -> None:
    """Plot the distribution and its approximation on a given axes."""
    axes.set_title(f'{property_report.params.name} probability plot')

    distribution_pair = property_report.distributions
    plot_probability_plot(distribution_pair, axes)


def _plot_qq_plot(property_report: ScalarNetworkPropertyReport, axes: plt.Axes) -> None:
    """Plot the distribution and its approximation on a given axes."""
    axes.set_title(f'{property_report.params.name} Q-Q plot')

    distribution_pair = property_report.distributions
    plot_qq_plot(distribution_pair, axes)
