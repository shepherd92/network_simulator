#!/usr/bin/env python3
"""Represent an empirical distribution based on sample values."""

from typing import Callable

import numpy as np
import numpy.typing as npt
from scipy.stats import kstest, gaussian_kde

from distribution.distribution import Distribution


class EmpiricalDistribution(Distribution):
    """Represents an empirical distribution."""

    def __init__(self, value_sequence: list[int | float]) -> None:
        """Calculate degree distribution of the network up to a certain maximum degree."""
        super().__init__()
        self._value_sequence = np.array(value_sequence)

        # filter out nan values
        self._value_sequence = self._value_sequence[~np.isnan(self._value_sequence)]

        if len(self._value_sequence) == 0:
            self._domain = Distribution.Domain(np.nan, np.nan)
            return

        self._value_sequence.sort()
        self._domain = Distribution.Domain(min(self._value_sequence), max(self._value_sequence))
        if len(np.unique(self._value_sequence)) != 1:
            self._density_estimator: Callable = gaussian_kde(self._value_sequence)
        else:
            self._density_estimator = self._dirac_density_estimator

        self._valid = True

    def get_value_sequence_in_domain(self, domain: Distribution.Domain) -> npt.NDArray[np.int_ | np.float_]:
        """Return the distribution values."""
        value_sequence_in_domain = self.value_sequence[
            (self.value_sequence >= domain.min_) &
            (self.value_sequence <= domain.max_)
        ]
        return value_sequence_in_domain

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Calculate the inverse CDF at the given values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'

        if self.empty:
            return np.full(quantiles_to_calculate.shape, np.nan)

        return np.quantile(self.value_sequence, quantiles_to_calculate)

    def test_normality(self) -> float:
        """Execute a standard normality test on the data and return the resulting p value."""
        return kstest(self.value_sequence, 'norm').pvalue

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the pdf values at x values provided that they are in the domain of the distribution."""
        if self.empty:
            return np.zeros_like(x_values)

        pdf_values = self._density_estimator(x_values)
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the cdf values at x values provided that they are in the domain of the distribution."""
        if self.empty:
            return np.full(x_values.shape, np.nan)
        cdf = np.ones_like(self.value_sequence).cumsum() / len(self.value_sequence)
        cdf_values = np.interp(x_values, self.value_sequence, cdf)

        return cdf_values

    def _dirac_density_estimator(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        result = np.zeros_like(x_values, dtype=np.float_)
        result[np.isclose(self.value_sequence[0], x_values)] = np.inf
        return result

    @property
    def natural_x_values(self) -> npt.NDArray[np.float_ | np.int_]:
        """Return some natural x values of the domain."""
        if self.empty:
            return np.empty((1))

        if self.all_values_are_the_same:
            return np.array([self.value_sequence[0]])

        if len(self.value_sequence < 100):
            return np.unique(self.value_sequence)

        num_of_bins = int(np.ceil(np.sqrt(len(self.value_sequence))))
        bin_edges = np.linspace(
            self.domain.min_,
            self.domain.max_,
            num_of_bins + 1,
            endpoint=True,
            dtype=float
        )
        bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
        bin_centers_in_domain = bin_centers[(bin_centers > self.domain.min_) & (bin_centers < self.domain.max_)]
        return bin_centers_in_domain

    @property
    def mean(self) -> np.float_:
        """Return the mean of the distribution."""
        return self.value_sequence.mean()

    @property
    def std_dev(self) -> np.float_:
        """Return the standard deviation of the distribution."""
        return self.value_sequence.std()

    @property
    def value_sequence(self) -> npt.NDArray[np.int_ | np.float_]:
        """Return the value sequence of the distribution."""
        return self._value_sequence

    @value_sequence.setter
    def value_sequence(self, value_sequence: npt.NDArray[np.int_ | np.float_]) -> None:
        self._value_sequence = value_sequence

    @property
    def empty(self) -> bool:
        """Return if the value sequence is empty."""
        return len(self.value_sequence) == 0

    @property
    def all_values_are_the_same(self) -> bool:
        """Tell if all values are the same in the value sequence."""
        if self.empty:
            return False
        return np.isclose(self.value_sequence, self.value_sequence[0]).all()

    def __str__(self) -> str:
        """Return string representation for reporting."""
        return f'Empirical Domain: {self.domain}'
