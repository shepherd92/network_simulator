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
    class DistributionParameters(TheoreticalDistribution.DistributionParameters):
        """Parameters of the uniform distribution."""

    @dataclass
    class DomainCalculation(TheoreticalDistribution.FittingParameters.DomainCalculation):
        """Parameters of domain calculation."""

    @dataclass
    class ParameterFitting(TheoreticalDistribution.FittingParameters.ParameterFitting):
        """Parameters of fitting."""

        class Method(Enum):
            """Method used for fitting the uniform distribution."""

            DEFAULT: int = auto()

        method: Method
        fixed_parameters: UniformDistribution.DistributionParameters

    def __init__(self) -> None:
        """Create a default uniform distribution."""
        super().__init__()
        self._parameters = UniformDistribution.DistributionParameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        return uniform.ppf(quantiles_to_calculate, *astuple(self.parameters))

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        result = super().info()
        result.update({
            'distribution_type': 'uniform',
        })

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_parameters: DomainCalculation,
    ) -> None:
        """Estimate the parameters of the uniform distribution."""
        self._domain.min_ = empirical_distribution.value_sequence.min()
        self._domain.max_ = empirical_distribution.value_sequence.max()

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        parameter_fitting_parameters: ParameterFitting
    ) -> None:
        """Estimate the parameters of the uniform distribution."""
        # value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        Method = UniformDistribution.ParameterFitting.Method
        if parameter_fitting_parameters.method == Method.DEFAULT:
            return
        else:
            assert False, f'Unknown fitting method: {parameter_fitting_parameters.method}.'

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        pdf_values = uniform.pdf(x_values, *astuple(self._parameters))
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = uniform.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @ property
    def parameters(self) -> UniformDistribution.DistributionParameters:
        """Return the parameters of the distribution."""
        return self._parameters
