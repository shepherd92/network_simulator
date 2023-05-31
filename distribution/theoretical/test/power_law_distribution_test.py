#!/usr/bin/env python3
"""Levy alpha-stable theoretical distribution test."""

import unittest

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pareto

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import create_fitting_parameters
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution

from reports.plotting_helper import plot_distribution_approximation


class PowerLawDistributionTest(unittest.TestCase):
    """Test case for testing the StableDistribution class."""

    def setUp(self):
        """Set up test case."""
        self.parameters = PowerLawDistribution.Parameters(exponent=2.5,)

        np.random.seed(seed=0)
        size = 100
        power_law_distributed_numbers = pareto.rvs(
            self.parameters.exponent - 1.,
            size=size,
        )
        assert power_law_distributed_numbers is not None
        empirical_distribution = EmpiricalDistribution(power_law_distributed_numbers.tolist())

        self.approximation = DistributionApproximation(empirical_distribution, TheoreticalDistribution.Type.POWER_LAW)
        fitting_parameters = create_fitting_parameters(TheoreticalDistribution.Type.POWER_LAW)
        self.approximation.fit(fitting_parameters)

    def test_fitting(self):
        """Test if the fitting method gives a reasonably good fit."""
        kolmogorov_smirnov_threshold = 0.1

        test_results = self.approximation.run_test()
        self._plot_pdfs()

        self.assertLess(
            test_results.kolmogorov_smirnov,
            kolmogorov_smirnov_threshold,
            f'Kolmogorov-Smirnov statistic is {test_results.kolmogorov_smirnov}, ' +
            f'which should be less than {kolmogorov_smirnov_threshold}.'
        )

    def _plot_pdfs(self):

        figure, axes = plt.subplots(1, 1)
        plot_distribution_approximation(self.approximation, axes)
        figure.show()


if __name__ == '__main__':

    unittest.main()
