#!/usr/bin/env python3
"""Levy alpha-stable theoretical distribution test."""

import unittest

import numpy as np
from scipy.stats import pareto

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import create_fitting_parameters
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class PowerLawDistributionTest(unittest.TestCase):
    """Test case for testing the StableDistribution class."""

    def setUp(self):
        """Set up test case."""
        self.parameters = PowerLawDistribution.Parameters(exponent=2.5,)

        np.random.seed(seed=0)
        size = 100000
        power_law_distributed_numbers = pareto.rvs(
            self.parameters.exponent - 1.,
            size=size,
        )
        assert power_law_distributed_numbers is not None
        self.empirical_distribution = EmpiricalDistribution(power_law_distributed_numbers.tolist())

    def test_fitting(self):
        """Test if the fitting method gives a reasonably good fit."""
        approximation = DistributionApproximation(self.empirical_distribution, TheoreticalDistribution.Type.POWER_LAW)
        fitting_parameters = create_fitting_parameters(TheoreticalDistribution.Type.POWER_LAW)
        fitting_parameters.fitting_method = PowerLawDistribution.FittingMethod.MAXIMUM_LIKELIHOOD_MLE_DOMAIN
        approximation.fit(fitting_parameters)

        self.assertAlmostEqual(approximation.theoretical.parameters.exponent, 2.5, delta=0.2)


if __name__ == '__main__':

    unittest.main()
