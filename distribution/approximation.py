#!/usr/bin/env python3
"""Theoretical approximation of an empirical distribution."""

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd

from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.factory import create_theoretical_distribution


@dataclass
class DistributionApproximation:
    """Encloses an empirical degree distribution and its approximation."""

    class TestResult(NamedTuple):
        """Results of the statistical tests for a scalar network property."""

        probability_plot_r_value: float = np.nan
        point_p_value: float = np.nan
        kolmogorov_smirnov: float = np.nan

        def __str__(self) -> str:
            """Return string representation for reporting."""
            result = f'probability_plot_r_value: {self.probability_plot_r_value:.4f}\n' + \
                f'Kolmogorov-Smirnov statistic: {self.kolmogorov_smirnov:.4f}'

            if not np.isnan(self.point_p_value):
                result += f'\npoint p-value: {self.point_p_value:.4f}'
            return result

    def __init__(self, empirical: EmpiricalDistribution, theoretical_type: TheoreticalDistribution.Type) -> None:
        """Construct an empirical theoretical distribution pair."""
        self._empirical = empirical
        self._theoretical = create_theoretical_distribution(theoretical_type)

    def fit(self):
        """Fit the theoretical distribution to the emirical distribution."""
        self.theoretical.fit(self.empirical)

    def run_test(self, point_value: float = np.nan) -> TestResult:
        """Run statistical tests."""
        if not self.valid:
            return DistributionApproximation.TestResult()

        kolmogorov_smirnov = self.theoretical.kolmogorov_smirnov(
            self.empirical,
            self.empirical.natural_x_values
        )
        probability_plot_r_value = self.calculate_probability_plot_r_value()

        if not np.isnan(point_value):
            p_value = self.theoretical.calc_p_values(np.array([point_value]))[0]
        else:
            p_value = np.nan

        test_result = DistributionApproximation.TestResult(
            probability_plot_r_value=probability_plot_r_value,
            point_p_value=p_value,
            kolmogorov_smirnov=kolmogorov_smirnov,
        )

        return test_result

    def generate_probability_plot_points(self) -> npt.NDArray[np.float_]:
        """Return the points of the probability points."""
        domain_intersection = self.theoretical.domain.intersect(self.empirical.domain)
        x_values = np.linspace(domain_intersection.min_, domain_intersection.max_, 100, endpoint=True)

        theoretical_cdf = self.theoretical.conditional_cdf(x_values)
        empirical_cdf = self.empirical.conditional_cdf(x_values)
        points = np.c_[empirical_cdf, theoretical_cdf]
        return points

    def generate_qq_plot_points(self) -> npt.NDArray[np.float_]:
        """Return the points of the probability points."""
        number_of_quantiles = len(self.empirical.natural_x_values)
        quantiles_to_calculate = np.linspace(0., 1., number_of_quantiles, endpoint=True)
        theoretical_quantiles = self.theoretical.calc_quantiles(quantiles_to_calculate)
        empirical_quantiles = self.empirical.calc_quantiles(quantiles_to_calculate)
        points = np.c_[theoretical_quantiles, empirical_quantiles]
        mean = self.empirical.value_sequence.mean()
        std = self.empirical.value_sequence.std()
        return (points - mean) / std

    def calculate_probability_plot_r_value(self) -> float:
        """Calculate the goodness of fit of a straight line fitted to the probability plot.

        The line to which we compare the points is y=x.
        """
        probability_plot_points = self.generate_probability_plot_points()
        ss_res = np.sum((probability_plot_points[:, 1] - probability_plot_points[:, 0])**2)
        probability_plot_r_value = 1 - ss_res
        return probability_plot_r_value

    def calculate_qq_plot_r_value(self) -> float:
        """Calculate the goodness of fit of a straight line fitted to the probability plot.

        The line to which we compare the points is y=x.
        """
        probability_plot_points = self.generate_qq_plot_points()
        ss_res = np.sum((probability_plot_points[:, 1] - probability_plot_points[:, 0])**2)
        ss_tot = np.sum((probability_plot_points[:, 1] - probability_plot_points[:, 1].mean())**2)
        probability_plot_r_value = 1 - ss_res / ss_tot
        return probability_plot_r_value

    def get_pdfs(self) -> pd.DataFrame:
        """Get the pdfs of the two distributions."""
        num_of_points = 1000
        x_values = np.linspace(
            self.empirical.domain.min_,
            self.empirical.domain.max_,
            num_of_points,
            endpoint=True
        )
        empirical_pdf = self.empirical.pdf(x_values)
        theoretical_pdf = self.theoretical.pdf(x_values)
        pdfs = pd.DataFrame(
            np.c_[empirical_pdf, theoretical_pdf],
            columns=['empirical', 'theoretical'],
            index=x_values,
        )
        return pdfs

    def get_confidence_intervals(self, confidence_levels: list[float]) -> pd.DataFrame:
        """Get confidence intervals for the specified confidence levels."""
        confidence_intervals = [
            self.empirical.calc_confidence_interval(confidence_level) +
            self.theoretical.calc_confidence_interval(confidence_level)
            for confidence_level in confidence_levels
        ]
        confidence_intervals_df = pd.DataFrame(
            confidence_intervals,
            columns=['empirical_lower', 'empirical_upper', 'theoretical_lower', 'theoretical_upper'],
            index=confidence_levels
        )
        return confidence_intervals_df

    @property
    def empirical(self) -> EmpiricalDistribution:
        """Return the empirical distribution part."""
        return self._empirical

    @property
    def theoretical(self) -> TheoreticalDistribution:
        """Return the theoretical approximation part."""
        return self._theoretical

    @property
    def valid(self) -> bool:
        """Return if both distributions are valid."""
        return self._theoretical.valid and self.empirical.valid

    def __str__(self) -> str:
        """Return a string of the main information for creating reports."""
        text = '\n'.join((
            str(self.theoretical),
            str(self.empirical)
        ))
        return text
