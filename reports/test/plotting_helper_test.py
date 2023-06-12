#!/usr/bin/env python3
"""Levy alpha-stable theoretical distribution test."""

import unittest

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.stats import levy_stable

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import create_default_fitting_parameters
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution

from reports.plotting_helper import (
    plot_distribution_approximation,
    plot_empirical_distribution_histogram_with_info,
)


class PlottingHelperTest(unittest.TestCase):
    """Test case for testing the StableDistribution class."""

    def setUp(self):
        """Set up test case."""
        np.random.seed(1)

    def test_plot_distribution_pdf_pair(self) -> None:
        """Test if plotting is correct."""
        figure, ((axes_poisson, axes_power_law), (axes_normal, axes_stable)) = plt.subplots(2, 2)
        figure.set_size_inches(20, 10)

        axes_poisson.set_title('Poisson Distribution')
        plot_distribution_approximation(
            PlottingHelperTest._get_dist_approximation(TheoreticalDistribution.Type.POISSON),
            axes_poisson
        )

        axes_power_law.set_title('Power law Distribution')
        plot_distribution_approximation(
            PlottingHelperTest._get_dist_approximation(TheoreticalDistribution.Type.POWER_LAW),
            axes_power_law
        )

        axes_normal.set_title('Normal Distribution')
        plot_distribution_approximation(
            PlottingHelperTest._get_dist_approximation(TheoreticalDistribution.Type.NORMAL),
            axes_normal
        )

        axes_stable.set_title('Stable Distribution')
        plot_distribution_approximation(
            PlottingHelperTest._get_dist_approximation(TheoreticalDistribution.Type.STABLE),
            axes_stable
        )

        figure.tight_layout()
        plt.show()
        plt.cla()

    def test_plot_empirical_pdf(self) -> None:
        """Test if plotting is correct."""
        figure, ((axes_poisson, axes_power_law), (axes_normal, axes_stable)) = plt.subplots(2, 2)
        figure.set_size_inches(20, 10)

        axes_poisson.set_title('Poisson Distribution')
        plot_empirical_distribution_histogram_with_info(
            PlottingHelperTest._get_empirical_dist(TheoreticalDistribution.Type.POISSON),
            EmpiricalDistribution.HistogramType.LINEAR,
            axes=axes_poisson
        )

        axes_normal.set_title('Normal Distribution')
        plot_empirical_distribution_histogram_with_info(
            PlottingHelperTest._get_empirical_dist(TheoreticalDistribution.Type.NORMAL),
            EmpiricalDistribution.HistogramType.LINEAR,
            axes=axes_normal
        )

        axes_stable.set_title('Stable Distribution')
        plot_empirical_distribution_histogram_with_info(
            PlottingHelperTest._get_empirical_dist(TheoreticalDistribution.Type.STABLE),
            EmpiricalDistribution.HistogramType.LINEAR,
            axes=axes_stable
        )

        axes_power_law.set_title('Power law Distribution')
        plot_empirical_distribution_histogram_with_info(
            PlottingHelperTest._get_empirical_dist(TheoreticalDistribution.Type.POWER_LAW),
            EmpiricalDistribution.HistogramType.LOGARITHMIC,
            axes=axes_power_law
        )

        figure.tight_layout()
        plt.show()
        plt.cla()

    def test_bin_sizes(self) -> None:
        """Test if plotting is correct."""
        figure, axes = plt.subplots(1, 1)

        plot_empirical_distribution_histogram_with_info(
            PlottingHelperTest._get_empirical_dist(TheoreticalDistribution.Type.UNIFORM),
            EmpiricalDistribution.HistogramType.LINEAR, axes
        )
        plt.show()
        plt.cla()

    @staticmethod
    def _get_dist_approximation(type_: TheoreticalDistribution.Type) -> DistributionApproximation:

        approximation = DistributionApproximation(PlottingHelperTest._get_empirical_dist(type_), type_)
        fitting_params = create_default_fitting_parameters(type_)
        approximation.fit(fitting_params)
        return approximation

    @ staticmethod
    def _get_empirical_dist(type_: TheoreticalDistribution.Type) -> EmpiricalDistribution:
        SAMPLE_SIZE = 1000

        def generate_power_law_in_domain(
            min_: float,
            max_: float,
            exponent: float,
            size: int = 1
        ) -> npt.NDArray[np.float_]:
            """Generate a power-law sample for pdf x^(g-1) in the interval [min_, max_]."""
            random = np.random.random(size=size)
            min_to_exponent, max_to_exponent = min_**exponent, max_**exponent
            return (min_to_exponent + (max_to_exponent - min_to_exponent) * random)**(1. / exponent)

        if type_ == TheoreticalDistribution.Type.UNIFORM:
            empirical_distribution = EmpiricalDistribution(np.random.uniform(1, 100, size=SAMPLE_SIZE).tolist())
        elif type_ == TheoreticalDistribution.Type.POISSON:
            empirical_distribution = EmpiricalDistribution(np.random.poisson(1, size=SAMPLE_SIZE).tolist())
        elif type_ == TheoreticalDistribution.Type.NORMAL:
            empirical_distribution = EmpiricalDistribution(
                np.random.normal(loc=3.0, scale=2.0, size=SAMPLE_SIZE).tolist()
            )
        elif type_ == TheoreticalDistribution.Type.STABLE:
            empirical_distribution = EmpiricalDistribution(
                levy_stable.rvs(alpha=1.5, beta=1.0, loc=3.0, scale=1.0, size=SAMPLE_SIZE).tolist()
            )
        elif type_ == TheoreticalDistribution.Type.POWER_LAW:
            empirical_distribution = EmpiricalDistribution(
                generate_power_law_in_domain(1, 1000, exponent=-1.5, size=SAMPLE_SIZE).tolist()
            )
        else:
            assert False, f'Unknown distribution type: {type_}'

        return empirical_distribution


if __name__ == '__main__':

    unittest.main()
