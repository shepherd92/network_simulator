#!/usr/bin/env python3
"""Power law theoretical distribution."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from logging import error

import numpy as np
import numpy.typing as npt

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class PowerLawDistribution(TheoreticalDistribution):
    """Power law theoretical distribution."""

    @dataclass
    class DistributionParameters(TheoreticalDistribution.DistributionParameters):
        """Parameters of the power law distribution."""

        exponent: float = np.nan

    @dataclass
    class DomainCalculation(TheoreticalDistribution.FittingParameters.DomainCalculation):
        """Parameters of domain calculation."""

        class Method(Enum):
            """Method used for the domain calculation."""

            LOGARITHMIC_RATIO = auto()
            MAXIMUM_LIKELIHOOD = auto()
            DETERMINISTIC = auto()
            QUANTILE = auto()

        method: Method

    @dataclass
    class MaximumLikelihoodDomain(DomainCalculation):
        """Maximum likelihood domain calculation parameters."""

        min_bounds: tuple[float, float] = (1., 100.)

    @dataclass
    class DeterministicDomain(DomainCalculation):
        """Deterministic domain calculation parameters."""

        min_: float = 10.
        max_: float = np.inf

    @dataclass
    class LogRatioDomain(DomainCalculation):
        """Logarithmic ratio domain calculation parameters."""

        min_ratio: float = 0.8
        max_ratio: float = 1.0

    @dataclass
    class QuantileDomain(DomainCalculation):
        """Quantile domain calculation parameters."""

        min_quantile: float = 0.8
        max_quantile: float = 1.0

    @dataclass
    class ParameterFitting(TheoreticalDistribution.FittingParameters.ParameterFitting):
        """Parameters of fitting."""

        class Method(Enum):
            """Method used for fitting the power law distribution."""

            LINEAR_REGRESSION = auto()
            MAXIMUM_LIKELIHOOD = auto()

        method: Method
        fixed_parameters: PowerLawDistribution.DistributionParameters

    def __init__(self) -> None:
        """Create a default power law distribution."""
        super().__init__()
        self._parameters = PowerLawDistribution.DistributionParameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        quantiles = self.domain.min_ * (1. - quantiles_to_calculate)**(1. / (1 - self._parameters.exponent))
        return quantiles

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        return {
            'distribution_type': 'power_law',
            'valid': self.valid,
            'domain': self.domain,
            'exponent': self._parameters.exponent,
        }

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_parameters: DomainCalculation,
    ) -> None:
        Method = PowerLawDistribution.DomainCalculation.Method
        if domain_calculation_parameters.method == Method.LOGARITHMIC_RATIO:
            self._determine_domain_logarithmic_ratio(
                empirical_distribution,
                domain_calculation_parameters
            )
        elif domain_calculation_parameters.method == Method.QUANTILE:
            self._determine_domain_quantile(
                empirical_distribution,
                domain_calculation_parameters
            )
        elif domain_calculation_parameters.method == Method.MAXIMUM_LIKELIHOOD:
            self._determine_domain_mle(
                empirical_distribution,
                domain_calculation_parameters
            )
        elif domain_calculation_parameters.method == Method.DETERMINISTIC:
            self._determine_domain_deterministic(
                empirical_distribution,
                domain_calculation_parameters
            )
        else:
            assert False, f'Unknown domain calculation method: {domain_calculation_parameters.method}.'

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        parameter_fitting_parameters: ParameterFitting,
    ) -> None:
        """Calculate the power law exponent of the distribution.

        Notations are consistent with the methods described in
        Power-Law Distributions in Empirical Data
        see: https://arxiv.org/pdf/0706.1062.pdf
        """
        if not np.isnan(parameter_fitting_parameters.fixed_parameters.exponent):
            self._parameters = parameter_fitting_parameters.fixed_parameters
            return

        Method = PowerLawDistribution.ParameterFitting.Method

        if parameter_fitting_parameters.method == Method.MAXIMUM_LIKELIHOOD:
            self._parameters.exponent = self._estimate_exponent_mle(
                empirical_distribution,
                parameter_fitting_parameters
            )
        elif parameter_fitting_parameters.method == Method.LINEAR_REGRESSION:
            self._parameters.exponent = self._estimate_exponent_linear_regression(
                empirical_distribution,
                parameter_fitting_parameters
            )
        else:
            assert False, f'Unknown fitting method: {parameter_fitting_parameters.method}.'

    def _determine_domain_deterministic(
        self,
        _: EmpiricalDistribution,
        domain_calculation_parameters: DeterministicDomain,
    ) -> None:
        """Find the domain at which the power-law holds."""
        self._domain.min_ = domain_calculation_parameters.min_
        self._domain.max_ = domain_calculation_parameters.max_

    def _determine_domain_mle(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_params: MaximumLikelihoodDomain,
    ) -> None:
        """Find the domain at which the power-law holds.

        min_bounds: minimum and maximum vlues of the domain minimum
        """
        min_bounds = Distribution.Domain(*domain_calculation_params.min_bounds)
        low = max(empirical_distribution.domain.min_, 1.)
        high = empirical_distribution.domain.max_
        valid_min_bounds = min_bounds.intersect(Distribution.Domain(low, high**0.75 * low**(1 - 0.75)))
        self._domain.max_ = empirical_distribution.domain.max_
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

    def _determine_domain_logarithmic_ratio(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_params: LogRatioDomain,
    ) -> None:
        """Find the minimum degree from which the power-law holds."""
        min_ratio = domain_calculation_params.min_ratio
        max_ratio = domain_calculation_params.max_ratio
        self._domain.min_ = empirical_distribution.domain.max_ ** min_ratio
        self._domain.max_ = empirical_distribution.domain.max_ ** max_ratio

    def _determine_domain_quantile(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_params: QuantileDomain,
    ) -> None:
        """Find the minimum degree from which the power-law holds."""
        min_quantile = domain_calculation_params.min_quantile
        max_quantile = domain_calculation_params.max_quantile
        self._domain.min_ = np.quantile(empirical_distribution.value_sequence, min_quantile)
        self._domain.max_ = np.quantile(empirical_distribution.value_sequence, max_quantile)

    def _estimate_exponent_mle(
        self,
        empirical_distribution: EmpiricalDistribution,
        _: ParameterFitting,
    ) -> DistributionParameters:
        """Estimate the power law  with maximum likelihood estimation.

        Notations are consistent with the methods described in
        Power-Law Distributions in Empirical Data
        see: https://arxiv.org/pdf/0706.1062.pdf
        """
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        if len(value_sequence) == 0:
            # warning(f'Value sequence is empty in domain [{self.domain.min_}, {self.domain.max_}].')
            return np.nan

        estimate = 1. + len(value_sequence) / (sum(np.log(value_sequence / (self.domain.min_ - 0.5))))
        return estimate

    def _estimate_exponent_linear_regression(
        self,
        empirical_distribution: EmpiricalDistribution,
        _: ParameterFitting,
    ) -> float:
        """Estimate the power law exponent."""
        value_counts = empirical_distribution.calc_value_counts()
        x_values = value_counts[:, 0]
        y_values = value_counts[:, 1]

        mask = \
            (x_values > 1e-10) & \
            (y_values > 1e-10) & \
            (x_values >= self.domain.min_) & \
            (x_values <= self.domain.max_)

        x_values = x_values[mask]
        y_values = y_values[mask]

        if len(x_values) < 3:
            error(f'Value counts is {len(x_values)} in domain [{self.domain.min_}, {self.domain.max_}].')
            return np.nan
        assert (x_values > 0.).all(), 'Non positive values occured in power law distribution.'

        coefficients = np.polyfit(np.log(x_values), np.log(y_values), 1)

        return - coefficients[0]

    def _objective_function(self, empirical_distribution: EmpiricalDistribution, guess: int) -> float:
        x_values = np.unique(empirical_distribution.value_sequence)
        if max(x_values) <= guess:
            return np.nan

        self._domain.min_ = guess
        self._parameters.exponent = self._estimate_exponent_mle(empirical_distribution)
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
        cdf_values = 1. - (x_values / self.domain.min_)**(1 - self.parameters.exponent)
        return cdf_values

    @property
    def parameters(self) -> PowerLawDistribution.DistributionParameters:
        """Return the parameters of the distribution."""
        return self._parameters
