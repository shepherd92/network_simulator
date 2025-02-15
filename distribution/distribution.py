#!/usr/bin/env python3
"""Represent the highest level base class for distributions."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd


class Distribution:
    """Represents a degree distribution."""

    @dataclass
    class Domain:
        """Domain of the degree distribution."""

        min_: float
        max_: float

        def intersect(self, other: Distribution.Domain) -> Distribution.Domain:
            """Return the intersection of two domains."""
            domain_min = max(self.min_, other.min_)
            domain_max = min(self.max_, other.max_)
            return Distribution.Domain(domain_min, domain_max)

        def filter_values(self, values: npt.NDArray[np.float64 | np.int_]) -> npt.NDArray[np.float64 | np.int_]:
            """Filter for those values that are in the domain."""
            values_in_domain = values[(values >= self.min_) & (values <= self.max_)]
            return values_in_domain

        @property
        def length(self) -> float:
            """Return the length of the interval."""
            return self.max_ - self.min_

        @property
        def valid(self) -> bool:
            """Return if the domain is valid."""
            nan_endpoints = np.isnan(self.min_) or np.isnan(self.max_)
            return not nan_endpoints and self.length >= 0

        def __str__(self) -> str:
            """Return string representation for reporting."""
            return f'[{float(self.min_):.4f}, {float(self.max_):.4f}]'

    def __init__(self) -> None:
        """Construct an invalid distribution."""
        self._domain = Distribution.Domain(np.nan, np.nan)
        self._valid: bool = False

    def pdf(self, x_values: npt.NDArray[np.float64 | np.int_]) -> npt.NDArray[np.float64]:
        """Return the PDF of the distribution evaluated at different x_values."""
        pdf_values = np.zeros_like(x_values, dtype=float)  # outside of the domain, the pdf is 0

        x_values_in_domain = x_values[(x_values >= self._domain.min_) & (x_values <= self._domain.max_)]
        pdf_values_in_domain = self._pdf_in_domain(x_values_in_domain)
        pdf_values[(x_values >= self._domain.min_) & (x_values <= self._domain.max_)] = pdf_values_in_domain

        pdf_values[np.isnan(x_values)] = np.nan

        return pdf_values

    def cdf(self, x_values: npt.NDArray[np.float64 | np.int_]) -> npt.NDArray[np.float64]:
        """Return the CDF of the distribution."""
        cdf_values = np.zeros_like(x_values, dtype=float)  # below domain min, the cdf is 0
        cdf_values[x_values > self._domain.max_] = 1.  # above domain min, the cdf is 1

        x_values_in_domain = x_values[(x_values >= self._domain.min_) & (x_values <= self._domain.max_)]
        cdf_values_in_domain = self._cdf_in_domain(x_values_in_domain)
        cdf_values[(x_values >= self._domain.min_) & (x_values <= self._domain.max_)] = cdf_values_in_domain

        cdf_values[np.isnan(x_values)] = np.nan
        return cdf_values

    def conditional_cdf(self, x_values: npt.NDArray[np.float64 | np.int_]) -> npt.NDArray[np.float64]:
        """Return the CDF conditioned on the variable is in the domain of x_values."""
        if not self.valid:
            return np.full(x_values.shape, np.nan)
        cdf = self.cdf(x_values)
        cdf_range = cdf.max() - cdf.min()
        conditional_cdf = (cdf - cdf.min()) / cdf_range
        return conditional_cdf

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Calculate the inverse CDF at the given values."""
        raise NotImplementedError

    def calc_confidence_interval(self, confidence: float) -> list[float]:
        """Calculcte the confidence interval of the distribution."""
        assert confidence < 1., f'The confidence level is {confidence} but it must be less than 1.'
        lower_quantile = 0.5 * (1. - confidence)
        upper_quantile = 1. - lower_quantile
        confidence_interval = self.calc_quantiles(np.array([lower_quantile, upper_quantile]))
        return confidence_interval.tolist()

    def kolmogorov_smirnov(self, other: Distribution, x_values: npt.NDArray[np.float64]) -> float:
        """Calculate the Kolmogorov-Smirnov statistic of two degree distributions."""
        if not self.valid or not other.valid:
            # only calculate the statistic if both distributions are valid
            return np.nan

        domain_intersection = self.domain.intersect(other.domain)
        x_values_in_domain = domain_intersection.filter_values(x_values)

        if len(x_values_in_domain) < 2:
            return np.nan

        this_cdf = self.conditional_cdf(x_values_in_domain)
        that_cdf = other.conditional_cdf(x_values_in_domain)

        kolmogorov_smirnov_statistic = np.max(np.abs(this_cdf - that_cdf))

        return kolmogorov_smirnov_statistic

    def calc_p_values(self, test_values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Calculate the two tailed p value for the given test value.

        Test if the test_value is from this distribution.
        """
        if not self.valid:
            return np.full(test_values.shape, np.nan)

        p_values = 2. * np.nanmin(np.c_[self.cdf(test_values), 1. - self.cdf(test_values)], axis=1)
        return p_values

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        return {
            'valid': self.valid,
            'domain_min': self.domain.min_,
            'domain_max': self.domain.max_,
        }

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        info = self.info()
        data_frame = pd.DataFrame(info, index=[0])
        data_frame.to_csv(save_path, index=False)

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float64 | np.int_]) -> npt.NDArray[np.float64]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        raise NotImplementedError

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float64 | np.int_]) -> npt.NDArray[np.float64]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        pdf_in_domain = self._pdf_in_domain(x_values)
        return pdf_in_domain.cumsum()

    @property
    def valid(self) -> bool:
        """Return if the degree distribution has valid values."""
        return self._valid

    @property
    def domain(self) -> Domain:
        """Return if the degree distribution has valid values."""
        return self._domain

    def __str__(self) -> str:
        """Return a string representation based on the network properties."""
        return '\n'.join([
            f'{key}: {item}'
            for key, item in self.info().items()
        ])
