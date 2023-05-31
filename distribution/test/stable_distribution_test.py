#!/usr/bin/env python3
"""Levy alpha-stable theoretical distribution test."""

import unittest

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import levy_stable

from distribution.approximation import DistributionApproximation
from distribution.factory import create_fitting_parameters
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.stable_distribution import StableDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class StableDistributionTest(unittest.TestCase):
    """Test case for testing the StableDistribution class."""

    def setUp(self):
        """Set up test case."""
        self.parameters = StableDistribution.Parameters(
            alpha=1.5,
            beta=1.,
            location=1000.,
            scale=100.,
        )

        np.random.seed(seed=0)

        size = 100

        stable_distributed_numbers = levy_stable.rvs(
            alpha=self.parameters.alpha,
            beta=self.parameters.beta,
            loc=self.parameters.location,
            scale=self.parameters.scale,
            size=size,
        )
        assert stable_distributed_numbers is not None
        self.empirical_distribution = EmpiricalDistribution(stable_distributed_numbers.tolist())

    def test_fitting(self):
        """Test if the fitting method gives a reasonably good fit."""
        kolmogorov_smirnov_threshold = 0.1

        approximation = DistributionApproximation(self.empirical_distribution, TheoreticalDistribution.Type.STABLE)
        fitting_parameters = create_fitting_parameters(TheoreticalDistribution.Type.STABLE)
        approximation.fit(fitting_parameters)

        test_results = approximation.run_test()
        # self._plot_pdfs()

        self.assertLess(
            test_results.kolmogorov_smirnov,
            kolmogorov_smirnov_threshold,
            f'Kolmogorov-Smirnov statistic is {test_results.kolmogorov_smirnov}, ' +
            f'which should be less than {kolmogorov_smirnov_threshold}.'
        )

    def _plot_pdfs(self):

        x_min = levy_stable.ppf(
            0.001,
            alpha=self.parameters.alpha,
            beta=self.parameters.beta,
            loc=self.parameters.location,
            scale=self.parameters.scale,
        )
        x_max = levy_stable.ppf(
            0.999,
            alpha=self.parameters.alpha,
            beta=self.parameters.beta,
            loc=self.parameters.location,
            scale=self.parameters.scale,
        )
        x_values = np.linspace(x_min, x_max, 1000, endpoint=True)

        empirical_pdf = self.approximation.empirical.pdf(x_values)
        theoretical_pdf = self.approximation.theoretical.pdf(x_values)
        plt.plot(x_values, empirical_pdf)
        plt.plot(x_values, theoretical_pdf)
        plt.show()


if __name__ == '__main__':

    unittest.main()
