#!/usr/bin/env python3
"""Poisson theoretical distribution."""

from __future__ import annotations

from dataclasses import dataclass, astuple
from enum import Enum, auto
from typing import NewType

import numpy as np
import numpy.typing as npt
from scipy.stats import poisson

from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


EmpiricalDistribution = NewType('EmpiricalDistribution', None)


class PoissonDistribution(TheoreticalDistribution):
    """Poisson theoretical degree distribution."""

    @dataclass
    class DistributionParameters(TheoreticalDistribution.DistributionParameters):
        """Parameters of the Poisson distribution."""

        lambda_: float = np.nan

    @dataclass
    class DomainCalculation(TheoreticalDistribution.FittingParameters.DomainCalculation):
        """Parameters of domain calculation."""

    @dataclass
    class ParameterFitting(TheoreticalDistribution.FittingParameters.ParameterFitting):
        """Parameters of fitting."""

        class Method(Enum):
            """Method used for fitting the Poisson distribution."""

            MAXIMUM_LIKELIHOOD: int = auto()

        method: Method
        fixed_parameters: PoissonDistribution.DistributionParameters

    def __init__(self) -> None:
        """Create a default Poisson distribution."""
        super().__init__()
        self._parameters = PoissonDistribution.DistributionParameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Calculate the quantiles."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quantiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return poisson.ppf(quantiles_to_calculate, mu=self.parameters.lambda_)

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        result = super().info()
        result.update({
            'distribution_type': 'poisson',
            'lambda': self.parameters.lambda_,
        })
        return result

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_parameters: DomainCalculation,
    ) -> None:
        self._domain.min_ = 0.
        self._domain.max_ = np.inf

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        parameter_fitting_parameters: ParameterFitting,
    ) -> None:
        """Calculate the parameter of the Poisson degree distribution."""
        if not np.isnan(parameter_fitting_parameters.fixed_parameters.lambda_):
            self._parameters = parameter_fitting_parameters.fixed_parameters
            return

        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)

        Method = PoissonDistribution.ParameterFitting.Method
        if parameter_fitting_parameters.method == Method.MAXIMUM_LIKELIHOOD:
            self._parameters = PoissonDistribution.DistributionParameters(value_sequence.mean())
        else:
            assert False, f'Unknown fitting method: {parameter_fitting_parameters.method}.'

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """PDF evaluated at the given x_values."""
        integers = np.array(list(range(int(x_values.min()), int(x_values.max()) + 1)))
        pmf_values = poisson.pmf(integers, *astuple(self._parameters))
        pdf_values = np.interp(x_values, integers, pmf_values)
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """CDF evaluated at the given x_values."""
        cdf_values = poisson.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @property
    def parameters(self) -> PoissonDistribution.DistributionParameters:
        """Return the parameters of the distribution."""
        return self._parameters
