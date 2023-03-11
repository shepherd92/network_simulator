#!/usr/bin/env python3
"""Poisson theoretical distribution."""

from __future__ import annotations

from dataclasses import dataclass, astuple

import numpy as np
import numpy.typing as npt
from scipy.stats import poisson

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class PoissonDistribution(TheoreticalDistribution):
    """Poisson theoretical degree distribution."""

    @dataclass
    class Parameters:
        """Parameters of the Poisson distribution."""

        lambda_: float = np.nan

    def __init__(self) -> None:
        """Create a default Poisson distribution."""
        super().__init__()
        self._parameters = PoissonDistribution.Parameters()
        self._domain = Distribution.Domain(0., np.inf)

    def _estimate_parameters(self, empirical_distribution: EmpiricalDistribution) -> None:
        """Calculate the parameter of the Poisson degree distribution."""
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        self._parameters = PoissonDistribution.Parameters(value_sequence.mean())

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return poisson.ppf(quantiles_to_calculate)

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        integers = np.array(list(range(int(x_values.min()), int(x_values.max()) + 1)))
        pmf_values = poisson.pmf(integers, *astuple(self._parameters))
        pdf_values = np.interp(x_values, integers, pmf_values)
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = poisson.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @property
    def parameters(self) -> PoissonDistribution.Parameters:
        """Return the parameters of the distribution."""
        return self._parameters

    def __str__(self) -> str:
        """Return string representation for reporting."""
        if not self.valid:
            return 'Invalid Theoretical Distribution'

        return '\n'.join([
            'Distribution: Poisson',
            f'Theoretical Domain: {self.domain}',
            f'Parameter: {self._parameters.lambda_:.4f}'
        ])
