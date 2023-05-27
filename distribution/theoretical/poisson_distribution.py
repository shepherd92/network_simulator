#!/usr/bin/env python3
"""Poisson theoretical distribution."""

from __future__ import annotations

from dataclasses import dataclass, astuple
from enum import Enum, auto

import numpy as np
import numpy.typing as npt
from scipy.stats import poisson

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class PoissonDistribution(TheoreticalDistribution):
    """Poisson theoretical degree distribution."""

    @dataclass
    class Parameters(TheoreticalDistribution.Parameters):
        """Parameters of the Poisson distribution."""

        lambda_: float = np.nan

    @dataclass
    class FittingParameters(TheoreticalDistribution.FittingParameters):
        """Parameters of how the fitting should be done."""

        fixed_parameters: PoissonDistribution.Parameters
        fitting_method: PoissonDistribution.FittingMethod

    class FittingMethod(Enum):
        """Method used for fitting the Poisson distribution."""

        MAXIMUM_LIKELIHOOD = auto()

    def __init__(self) -> None:
        """Create a default Poisson distribution."""
        super().__init__()
        self._parameters = PoissonDistribution.Parameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return poisson.ppf(quantiles_to_calculate)

    def get_info_as_dict(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        return {
            'distribution_type': 'poisson',
            'valid': self.valid,
            'domain_min': self.domain.min_,
            'domain_max': self.domain.max_,
            'lambda': self.parameters.lambda_,
        }

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters
    ) -> None:
        if fitting_parameters.fitting_method == PoissonDistribution.FittingMethod.MAXIMUM_LIKELIHOOD:
            self._domain = Distribution.Domain(0., np.inf)
        else:
            assert False, f'Unknown fitting method: {fitting_parameters.fitting_method}.'

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters
    ) -> None:
        """Calculate the parameter of the Poisson degree distribution."""
        if not np.isnan(fitting_parameters.fixed_parameters.lambda_):
            self._parameters = fitting_parameters
            return

        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)

        if fitting_parameters.fitting_method == PoissonDistribution.FittingMethod.MAXIMUM_LIKELIHOOD:
            self._parameters = PoissonDistribution.Parameters(value_sequence.mean())
        else:
            assert False, f'Unknown fitting method: {fitting_parameters.fitting_method}.'

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
