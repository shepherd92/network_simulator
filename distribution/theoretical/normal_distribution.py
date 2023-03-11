#!/usr/bin/env python3
"""Normal theoretical distribution."""

from __future__ import annotations

from dataclasses import astuple, dataclass
from logging import warning

import numpy as np
import numpy.typing as npt
from scipy.stats import norm

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class NormalDistribution(TheoreticalDistribution):
    """Normal theoretical distribution."""

    @dataclass
    class Parameters:
        """Parameters of the normal distribution."""

        mean: float = np.nan
        std: float = np.nan

    def __init__(self) -> None:
        """Create a default normal distribution."""
        super().__init__()
        self._parameters = NormalDistribution.Parameters()
        self._domain = Distribution.Domain(-np.inf, np.inf)

    def _estimate_parameters(self, empirical_distribution: EmpiricalDistribution) -> None:
        """Estimate the parameters of the normal distribution."""
        self._parameters.mean = empirical_distribution.value_sequence.mean()
        self._parameters.std = empirical_distribution.value_sequence.std()

        if np.isclose(self._parameters.std, 0.):
            warning('Normal distribution fitting is invalid as data has 0 standard deviation.')

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return norm.ppf(quantiles_to_calculate, *astuple(self.parameters))

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        pdf_values = norm.pdf(x_values, *astuple(self._parameters))
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = norm.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @property
    def parameters(self) -> NormalDistribution.Parameters:
        """Return the parameters of the distribution."""
        return self._parameters

    def __str__(self) -> str:
        """Return string representation for reporting."""
        if not self.valid:
            return 'Invalid Theoretical Distribution'

        return '\n'.join([
            'Distribution: Normal',
            f'Theoretical Domain: {self.domain}',
            f'Mean: {self._parameters.mean:.4f}',
            f'Std: {self._parameters.std:.4f}'
        ])
