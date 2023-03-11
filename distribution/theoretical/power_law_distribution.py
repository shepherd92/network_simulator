#!/usr/bin/env python3
"""Power law theoretical distribution."""

from __future__ import annotations

from dataclasses import dataclass
from logging import warning

import numpy as np
import numpy.typing as npt

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class PowerLawDistribution(TheoreticalDistribution):
    """Power law theoretical distribution."""

    @dataclass
    class Parameters:
        """Parameters of the power law distribution."""

        exponent: float = np.nan

    def __init__(self) -> None:
        """Create a default power law distribution."""
        super().__init__()
        self._domain = Distribution.Domain(1., np.inf)
        self._parameters = PowerLawDistribution.Parameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        quantiles = self.domain.min_ * (1. - quantiles_to_calculate)**(1. / (1 - self._parameters.exponent))
        return quantiles

    def _estimate_parameters(self, empirical_distribution: EmpiricalDistribution) -> None:
        """Calculate the power law exponent of the distribution.

        Notations are consistent with the methods described in
        Power-Law Distributions in Empirical Data
        see: https://arxiv.org/pdf/0706.1062.pdf
        """
        self._determine_domain(empirical_distribution)
        self._parameters = self._estimate_exponent(empirical_distribution)

    def _determine_domain(self, empirical_distribution: EmpiricalDistribution) -> None:
        """Find the domain at which the power-law holds.

        min_bounds: minimum and maximum vlues of the domain minimum
        """
        # calculate the intersection of the domains so that at least 25% of the domain remains above the minimum
        min_bounds = Distribution.Domain(1, 100)
        low = max(empirical_distribution.domain.min_, 1.)
        high = empirical_distribution.domain.max_
        valid_min_bounds = min_bounds.intersect(Distribution.Domain(low, high**0.75 * low**(1 - 0.75)))
        # create an array with guesses in the 0th column, and respective Kolmogorov-Smirnov
        # statistics in the 1st column
        ks_statistics: npt.NDArray[np.float_] = np.array([
            [guess, self._objective_function(empirical_distribution, guess)]
            for guess in range(int(valid_min_bounds.min_), int(valid_min_bounds.max_) + 1)
        ])

        if len(ks_statistics) == 0 or np.isnan(ks_statistics[:, 1]).all():
            self._domain.min_ = np.nan
            return

        self._domain.min_ = np.round(ks_statistics[np.nanargmin([ks_statistics[:, 1]]), 0])

    def _estimate_exponent(self, empirical_distribution: EmpiricalDistribution) -> PowerLawDistribution.Parameters:
        """Estimate the power law exponent."""
        if not self.domain.valid:
            return PowerLawDistribution.Parameters()

        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        if len(value_sequence) == 0:
            warning(f'Value sequence is empty in domain [{self.domain.min_}, {self.domain.max_}].')
            return PowerLawDistribution.Parameters(np.nan)

        estimate = 1. + len(value_sequence) / (sum(np.log(value_sequence / (self.domain.min_ - 0.5))))
        return PowerLawDistribution.Parameters(estimate)

    def _objective_function(self, empirical_distribution: EmpiricalDistribution, guess: int) -> float:
        x_values = empirical_distribution.natural_x_values
        if max(x_values) <= guess:
            return np.nan

        self._domain.min_ = guess
        self._parameters = self._estimate_exponent(empirical_distribution)
        self._valid = True  # pylint: disable=attribute-defined-outside-init

        ks_statistic = self.kolmogorov_smirnov(empirical_distribution, x_values)
        return ks_statistic

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        x_min = self.domain.min_
        exponent = self._parameters.exponent
        pdf_values = (exponent - 1.) / x_min * (x_values / x_min)**(-exponent)
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = 1. - (x_values / self.domain.min_)**(1 - self._parameters.exponent)
        return cdf_values

    @property
    def parameters(self) -> PowerLawDistribution.Parameters:
        """Return the parameters of the distribution."""
        return self._parameters

    def __str__(self) -> str:
        """Return string representation for reporting."""
        if not self.valid:
            return 'Invalid Theoretical Distribution'

        return '\n'.join([
            'Distribution: Power Law',
            f'Theoretical Domain: {self.domain}',
            f'Exponent: {self._parameters.exponent:.4f}'
        ])
