#!/usr/bin/env python3
"""Base class for theoretical distributions."""

from __future__ import annotations

from dataclasses import dataclass, astuple
from enum import Enum, auto
from logging import warning

import numpy as np

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution


class TheoreticalDistribution(Distribution):
    """Represents a theoretical degree distribution."""

    class Type(Enum):
        """Type of the theoretical distribution."""

        UNIFORM: int = auto()
        POISSON: int = auto()
        POWER_LAW: int = auto()
        NORMAL: int = auto()
        STABLE: int = auto()
        INVALID: int = auto()

    @dataclass
    class DistributionParameters:
        """Represent the parameters of the distribution."""

    @dataclass
    class FittingParameters:
        """Parameters of how the fitting should be done."""

        @dataclass
        class DomainCalculation:
            """Parameters of domain calculation."""

        @dataclass
        class ParameterFitting:
            """Parameters of domain calculation."""

        domain_calculation: DomainCalculation
        parameter_fitting: ParameterFitting

    def __init__(self) -> None:
        """Create a default theoretical distribution."""
        super().__init__()
        self._parameters = TheoreticalDistribution.DistributionParameters()

    def fit(self, empirical_distribution: EmpiricalDistribution, fitting_parameters: FittingParameters) -> None:
        """Fit the parameters of the probability distribution."""
        if len(empirical_distribution.value_sequence) < 2:
            warning(
                f'Empirical degree distribution contains {len(empirical_distribution.value_sequence)} in the domain.')
            self._valid = False  # pylint: disable=attribute-defined-outside-init
            return

        self._fit_domain(empirical_distribution, fitting_parameters.domain_calculation)
        if not self.domain.valid:
            return

        self._fit_parameters(empirical_distribution, fitting_parameters.parameter_fitting)
        self._valid = not any(  # pylint: disable=attribute-defined-outside-init
            np.isnan(x)
            for x in astuple(self._parameters)
        ) and self.domain.valid

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        return {
            'valid': self.valid,
            'domain_min': self.domain.min_,
            'domain_max': self.domain.max_,
        }

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_parameters: FittingParameters.DomainCalculation,
    ) -> None:
        raise NotImplementedError

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        parameter_fitting_parameters: FittingParameters.ParameterFitting,
    ) -> None:
        raise NotImplementedError

    @property
    def parameters(self) -> DistributionParameters:
        """Return the parameters of the distribution."""
        return self._parameters

    @parameters.setter
    def parameters(self, parameters: DistributionParameters) -> None:
        self._parameters = parameters
