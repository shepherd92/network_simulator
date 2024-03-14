#!/usr/bin/env python3
"""Normal theoretical distribution."""

from __future__ import annotations

from dataclasses import astuple, dataclass
from enum import Enum, auto
from logging import warning

import numpy as np
import numpy.typing as npt
from scipy.stats import norm

from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class NormalDistribution(TheoreticalDistribution):
    """Normal theoretical distribution."""

    @dataclass
    class DistributionParameters(TheoreticalDistribution.DistributionParameters):
        """Parameters of the normal distribution."""

        mean: float = np.nan
        std: float = np.nan

    @dataclass
    class DomainCalculation(TheoreticalDistribution.FittingParameters.DomainCalculation):
        """Parameters of domain calculation."""

    @dataclass
    class ParameterFitting(TheoreticalDistribution.FittingParameters.ParameterFitting):
        """Parameters of fitting."""

        class Method(Enum):
            """Method used for fitting the normal distribution."""

            MAXIMUM_LIKELIHOOD: int = auto()
            MATCH_QUANTILE = auto()

        method: Method
        fixed_parameters: NormalDistribution.DistributionParameters

    @dataclass
    class ParameterFittingMatchQuantile(ParameterFitting):
        """Parameters of the match percentile fitting."""

        quantile: float

    def __init__(self) -> None:
        """Create a default normal distribution."""
        super().__init__()
        self._parameters = NormalDistribution.DistributionParameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return norm.ppf(quantiles_to_calculate, *astuple(self.parameters))

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        return {
            'distribution_type': 'normal',
            'valid': self.valid,
            'domain': self.domain,
            'mean': self.parameters.mean,
            'std': self.parameters.std,
        }

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_parameters: DomainCalculation,
    ) -> None:
        self._domain.min_ = -np.inf
        self._domain.max_ = +np.inf

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        parameter_fitting_parameters: ParameterFitting,
    ) -> None:
        """Estimate the parameters of the normal distribution."""
        Method = NormalDistribution.ParameterFitting.Method
        if parameter_fitting_parameters.method == Method.MAXIMUM_LIKELIHOOD:
            self._parameters = self._fit_maximum_likelihood(empirical_distribution, parameter_fitting_parameters)
        elif parameter_fitting_parameters.method == Method.MATCH_QUANTILE:
            self._parameters = self._fit_match_quantile(empirical_distribution, parameter_fitting_parameters)
        else:
            assert False, f'Unknown fitting method: {parameter_fitting_parameters.method}.'

        if np.isclose(self._parameters.std, 0.):
            warning('Normal distribution fitting is invalid as data has 0 standard deviation.')

    def _fit_maximum_likelihood(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: ParameterFitting
    ) -> DistributionParameters:
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        mean = value_sequence.mean() \
            if np.isnan(fitting_parameters.fixed_parameters.mean) else fitting_parameters.fixed_parameters.mean
        std = value_sequence.std() \
            if np.isnan(fitting_parameters.fixed_parameters.std) else fitting_parameters.fixed_parameters.std
        return NormalDistribution.DistributionParameters(mean, std)

    def _fit_match_quantile(
        self,
        empirical_distribution: EmpiricalDistribution,
        fitting_parameters: ParameterFittingMatchQuantile,
    ) -> DistributionParameters:
        """Set the variance so that the specified quantiles match.

        If F_data(x) = quantile_to_compute, then F_approximation(x) should be quantile_to_compute as well.
        F_A (F_D^-1 (q)) = q, but as F_A is the CDF of a normal distribution, we can standardize:
        q = F_Z( (F_D^-1 (q) - mu) / std ), so by expressing std:
        std = (F_D^-1 (q) - mu) / F_Z^-1 (q)
        """
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        mean = value_sequence.mean() \
            if np.isnan(fitting_parameters.fixed_parameters.mean) else fitting_parameters.fixed_parameters.mean

        quantile = fitting_parameters.quantile
        if np.isnan(fitting_parameters.fixed_parameters.std):
            quantile_of_data = empirical_distribution.calc_quantiles(np.array([quantile]))[0]
            std = (quantile_of_data - mean) / norm.ppf(quantile)
        else:
            std = fitting_parameters.fixed_parameters.std

        return NormalDistribution.DistributionParameters(mean, std)

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        pdf_values = norm.pdf(x_values, *astuple(self._parameters))
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = norm.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @property
    def parameters(self) -> NormalDistribution.DistributionParameters:
        """Return the parameters of the distribution."""
        return self._parameters
