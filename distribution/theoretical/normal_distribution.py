#!/usr/bin/env python3
"""Normal theoretical distribution."""

from __future__ import annotations

from dataclasses import astuple, dataclass
from enum import Enum, auto
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
    class Parameters(TheoreticalDistribution.Parameters):
        """Parameters of the normal distribution."""

        mean: float = np.nan
        std: float = np.nan

    @dataclass
    class FittingParameters(TheoreticalDistribution.FittingParameters):
        """Parameters of how the fitting should be done."""

        fixed_parameters: NormalDistribution.Parameters
        fitting_method: NormalDistribution.FittingMethod

    class FittingMethod(Enum):
        """Method used for fitting the normal distribution."""

        MAXIMUM_LIKELIHOOD = auto()
        MATCH_2_DOT_5_PERCENTILE = auto()

    def __init__(self) -> None:
        """Create a default normal distribution."""
        super().__init__()
        self._parameters = NormalDistribution.Parameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return norm.ppf(quantiles_to_calculate, *astuple(self.parameters))

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters
    ) -> None:
        if fitting_parameters.fitting_method in [
            NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
            NormalDistribution.FittingMethod.MATCH_2_DOT_5_PERCENTILE,
        ]:
            self._domain = Distribution.Domain(-np.inf, np.inf)
        else:
            assert False, f'Unknown fitting method: {fitting_parameters.fitting_method}.'

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters
    ) -> None:
        """Estimate the parameters of the normal distribution."""
        if fitting_parameters.fitting_method == NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD:
            self._parameters = self._fit_maximum_likelihood(empirical_distribution, fitting_parameters)
        elif fitting_parameters.fitting_method == NormalDistribution.FittingMethod.MATCH_2_DOT_5_PERCENTILE:
            self._parameters = self._fit_match_quantile(empirical_distribution, fitting_parameters, 0.025)
        else:
            assert False, f'Unknown fitting method: {fitting_parameters.fitting_method}.'

        if np.isclose(self._parameters.std, 0.):
            warning('Normal distribution fitting is invalid as data has 0 standard deviation.')

    def _fit_maximum_likelihood(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters
    ) -> Parameters:
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        mean = value_sequence.mean() \
            if np.isnan(fitting_parameters.fixed_parameters.mean) else fitting_parameters.fixed_parameters.mean
        std = value_sequence.std() \
            if np.isnan(fitting_parameters.fixed_parameters.std) else fitting_parameters.fixed_parameters.std
        return NormalDistribution.Parameters(mean, std)

    def _fit_match_quantile(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: FittingParameters,
        quantile_to_compute: float
    ) -> Parameters:
        """Set the variance so that the specified quantiles match.

        If F_data(x) = quantile_to_compute, then F_approximation(x) should be quantile_to_compute as well.
        F_A (F_D^-1 (q)) = q, but as F_A is the CDF of a normal distribution, we can standardize:
        q = F_Z( (F_D^-1 (q) - mu) / std ), so by expressing std:
        std = (F_D^-1 (q) - mu) / F_Z^-1 (q)
        """
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        mean = value_sequence.mean() \
            if np.isnan(fitting_parameters.fixed_parameters.mean) else fitting_parameters.fixed_parameters.mean

        if np.isnan(fitting_parameters.fixed_parameters.std):
            quantile_of_data = empirical_distribution.calc_quantiles(np.array([quantile_to_compute]))[0]
            std = (quantile_of_data - mean) / norm.ppf(quantile_to_compute)
        else:
            std = fitting_parameters.fixed_parameters.std

        return NormalDistribution.Parameters(mean, std)

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
