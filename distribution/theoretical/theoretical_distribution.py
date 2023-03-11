#!/usr/bin/env python3
"""Base class for theoretical distributions."""

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

        POISSON: int = auto()
        POWER_LAW: int = auto()
        NORMAL: int = auto()
        STABLE: int = auto()
        INVALID: int = auto()

    @dataclass
    class Parameters:
        """Represent the parameters of the distribution."""

    def __init__(self) -> None:
        """Create a default theoretical distribution."""
        super().__init__()
        self._parameters = TheoreticalDistribution.Parameters()

    def fit(self, empirical_distribution: EmpiricalDistribution) -> None:
        """Fit the parameters of the probability distribution."""
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)

        if len(empirical_distribution.value_sequence) < 2:
            warning(f'Empirical degree distribution contains {len(value_sequence)} in the domain.')
            self._valid = False  # pylint: disable=attribute-defined-outside-init
            return

        self._estimate_parameters(empirical_distribution)
        self._valid = not any(  # pylint: disable=attribute-defined-outside-init
            np.isnan(x)
            for x in astuple(self._parameters)
        ) and self.domain.valid

    def _estimate_parameters(self, empirical_distribution: EmpiricalDistribution) -> None:
        raise NotImplementedError

    @property
    def parameters(self) -> Parameters:
        """Return the parameters of the distribution."""
        return self._parameters
