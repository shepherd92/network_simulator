#!/usr/bin/env python3
"""Uniform theoretical distribution."""

from __future__ import annotations

from dataclasses import astuple, dataclass
from enum import Enum, auto

import numpy as np
import numpy.typing as npt
from scipy.stats import uniform

from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class UniformDistribution(TheoreticalDistribution):
    """Uniform theoretical distribution."""

    @dataclass
    class Parameters(TheoreticalDistribution.Parameters):
        """Parameters of the uniform distribution."""

    @dataclass
    class FittingParameters(TheoreticalDistribution.FittingParameters):
        """Parameters of how the fitting should be done."""

        fixed_parameters: UniformDistribution.Parameters
        fitting_method: UniformDistribution.FittingMethod

    class FittingMethod(Enum):
        """Method used for fitting the uniform distribution."""

        DEFAULT: int = auto()

    def __init__(self) -> None:
        """Create a default uniform distribution."""
        super().__init__()
        self._parameters = UniformDistribution.Parameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return uniform.ppf(quantiles_to_calculate, *astuple(self.parameters))

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters
    ) -> None:
        """Estimate the parameters of the uniform distribution."""
        if fitting_parameters.fitting_method == UniformDistribution.FittingMethod.DEFAULT:
            self._domain.min_ = empirical_distribution.value_sequence.min()
            self._domain.max_ = empirical_distribution.value_sequence.max()
        else:
            assert False, f'Unknown fitting method: {fitting_parameters.fitting_method}.'

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters
    ) -> None:
        """Estimate the parameters of the uniform distribution."""
        # value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        if fitting_parameters.fitting_method == UniformDistribution.FittingMethod.DEFAULT:
            return
        else:
            assert False, f'Unknown fitting method: {fitting_parameters.fitting_method}.'

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        pdf_values = uniform.pdf(x_values, *astuple(self._parameters))
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = uniform.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @property
    def parameters(self) -> UniformDistribution.Parameters:
        """Return the parameters of the distribution."""
        return self._parameters

    def __str__(self) -> str:
        """Return string representation for reporting."""
        if not self.valid:
            return 'Invalid Theoretical Distribution'

        return '\n'.join([
            'Distribution: Normal',
            f'Theoretical Domain: {self.domain}',
            f'Min: {self._parameters.min_:.4f}',
            f'Max: {self._parameters.max_:.4f}'
        ])
