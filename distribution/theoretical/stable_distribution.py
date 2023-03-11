#!/usr/bin/env python3
"""Levy alpha-stable theoretical distribution."""

from __future__ import annotations

from dataclasses import dataclass, astuple

import numpy as np
import numpy.typing as npt
from scipy.stats import levy_stable

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class StableDistribution(TheoreticalDistribution):
    """Levy-stable theoretical distribution."""

    @dataclass
    class Parameters(TheoreticalDistribution.Parameters):
        """Parameters of the stable distribution."""

        alpha: float = np.nan
        beta: float = np.nan
        location: float = np.nan
        scale: float = np.nan

    def __init__(self) -> None:
        """Create a default stable distribution."""
        super().__init__()
        self._domain = Distribution.Domain(-np.inf, np.inf)
        self._parameters = StableDistribution.Parameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        quantiles = levy_stable.ppf(quantiles_to_calculate, *astuple(self._parameters))
        return quantiles

    def _estimate_parameters(self, empirical_distribution: EmpiricalDistribution) -> None:
        def pconv(alpha: float, beta: float, location: float, sigma: float) -> StableDistribution.Parameters:
            result = StableDistribution.Parameters(
                alpha=alpha,
                beta=beta,
                location=location - sigma * beta * np.tan(np.pi * alpha / 2.),
                scale=sigma
            )
            return result

        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        self._parameters = pconv(*levy_stable._fitstart(value_sequence))  # pylint: disable=protected-access

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        pdf_values = levy_stable.pdf(x_values, *astuple(self._parameters))
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = levy_stable.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @property
    def parameters(self) -> StableDistribution.Parameters:
        """Return the parameters of the distribution."""
        return self._parameters

    def __str__(self) -> str:
        """Return string representation for reporting."""
        if not self.valid:
            return 'Invalid Theoretical Distribution'

        return '\n'.join([
            'Distribution: Stable',
            f'Theoretical Domain: {self.domain}',
            f'Alpha: {self._parameters.alpha:.4f}',
            f'Beta: {self._parameters.beta}',
            f'Location: {self._parameters.location}',
            f'Scale: {self._parameters.scale}'
        ])
