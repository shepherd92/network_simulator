#!/usr/bin/env python3
"""Theoretical approximation of an empirical distribution."""

from pathlib import Path
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd

from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.factory import create_theoretical_distribution


class DistributionApproximation:
    """Encloses an empirical degree distribution and its approximation."""

    class TestResult(NamedTuple):
        """Results of the statistical tests for a scalar network property."""

        probability_plot_r_value: float = np.nan
        qq_plot_r_value: float = np.nan
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
        self._type = theoretical_type
        self._empirical = empirical
        self._theoretical = create_theoretical_distribution(theoretical_type)
        self._pdfs: pd.DataFrame | None = None

    def fit(self, fitting_parameters: TheoreticalDistribution.FittingParameters) -> None:
        """Fit the theoretical distribution to the emirical distribution."""
        self.theoretical.fit(self.empirical, fitting_parameters)

    def run_test(self, point_value: float = np.nan) -> TestResult:
        """Run statistical tests."""
        if not self.valid:
            return DistributionApproximation.TestResult()

        kolmogorov_smirnov = self.theoretical.kolmogorov_smirnov(
            self.empirical,
            np.unique(self.empirical.value_sequence)
        )
        probability_plot_r_value = self.calculate_probability_plot_r_value()
        qq_plot_r_value = self.calculate_qq_plot_r_value()

        if not np.isnan(point_value):
            p_value = self.theoretical.calc_p_values(np.array([point_value]))[0]
        else:
            p_value = np.nan

        test_result = DistributionApproximation.TestResult(
            probability_plot_r_value=probability_plot_r_value,
            qq_plot_r_value=qq_plot_r_value,
            point_p_value=p_value,
            kolmogorov_smirnov=kolmogorov_smirnov,
        )

        return test_result

    def calculate_probability_plot_r_value(self) -> float:
        """Calculate the goodness of fit of a straight line fitted to the probability plot.

        The line to which we compare the points is y=x.
        """
        points = self.generate_probability_plot_points()
        return self._calc_r_squared_value(points)

    def generate_probability_plot_points(self) -> npt.NDArray[np.float_]:
        """Return the points of the probability points."""
        domain_intersection = self.theoretical.domain.intersect(self.empirical.domain)
        x_values = np.linspace(domain_intersection.min_, domain_intersection.max_, 100, endpoint=True)

        theoretical_cdf = self.theoretical.conditional_cdf(x_values)
        empirical_cdf = self.empirical.conditional_cdf(x_values)
        points = np.c_[empirical_cdf, theoretical_cdf]
        return points

    def calculate_qq_plot_r_value(self) -> float:
        """Calculate the goodness of fit of a straight line fitted to the probability plot.

        The line to which we compare the points is y=x.
        """
        points = self.generate_qq_plot_points()
        return self._calc_r_squared_value(points)

    def generate_qq_plot_points(self) -> npt.NDArray[np.float_]:
        """Return the points of the probability points."""
        number_of_quantiles = len(np.unique(self.empirical.value_sequence))
        quantiles_to_calculate = np.linspace(0., 1., number_of_quantiles, endpoint=True)
        theoretical_quantiles = self.theoretical.calc_quantiles(quantiles_to_calculate)
        empirical_quantiles = self.empirical.calc_quantiles(quantiles_to_calculate)
        points = np.c_[theoretical_quantiles, empirical_quantiles]
        mean = self.empirical.value_sequence.mean()
        std = self.empirical.value_sequence.std()
        return (points - mean) / std

    def save(self, save_directory: Path) -> None:
        """Save all information of the aroximation."""
        save_directory.mkdir(parents=True, exist_ok=True)
        self.save_info(save_directory / 'distribution_info.csv')

        if len(self.empirical.value_sequence) == 0:
            return

        self.pdfs.to_csv(save_directory / 'pdfs.csv', float_format='%.9f')
        np.savetxt(save_directory / 'value_sequence.csv', self.empirical.value_sequence, delimiter=',')
        if self.empirical.value_sequence.dtype == np.int_:
            value_counts = self.empirical.calc_value_counts()
            np.savetxt(save_directory / 'value_counts.csv', value_counts, delimiter=',')

        confidence_levels = [0.9, 0.95, 0.99]
        confidence_intervals = self.get_confidence_intervals(confidence_levels)
        confidence_intervals.to_csv(save_directory / 'confidence_intervals.csv', float_format='%.4f')

        probability_plot_points = self.generate_probability_plot_points()
        np.savetxt(save_directory / 'probability_plot_points.csv', probability_plot_points, delimiter=',')
        qq_plot_points = self.generate_qq_plot_points()
        np.savetxt(save_directory / 'qq_plot_points.csv', qq_plot_points, delimiter=',')

    def _get_pdfs(self) -> pd.DataFrame:
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

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        info = self.get_info_as_dict()
        data_frame = pd.DataFrame(info, index=[0])
        data_frame.to_csv(save_path, index=False)

    def get_info_as_dict(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        empirical_info = self.empirical.get_info_as_dict()
        theoretical_info = self.theoretical.get_info_as_dict()

        joint_info: dict[str, int | float] = {}
        for key, value in empirical_info.items():
            joint_info[f'empirical_{key}'] = value
        for key, value in theoretical_info.items():
            joint_info[f'theoretical_{key}'] = value

        test_result = self.run_test()
        joint_info['kolmogorov_smirnov'] = test_result.kolmogorov_smirnov
        joint_info['probability_plot_r_value'] = test_result.probability_plot_r_value
        joint_info['qq_plot_r_value'] = test_result.qq_plot_r_value

        return joint_info

    @staticmethod
    def _calc_r_squared_value(points: npt.NDArray) -> float:

        ss_res = np.sum((points[:, 1] - points[:, 0])**2)
        ss_tot = np.sum((points[:, 1] - points[:, 1].mean())**2)
        probability_plot_r_value = 1 - ss_res / ss_tot
        return probability_plot_r_value

    @property
    def type(self) -> TheoreticalDistribution.Type:
        """Return the theoretical approximation type."""
        return self._type

    @property
    def empirical(self) -> EmpiricalDistribution:
        """Return the empirical distribution part."""
        return self._empirical

    @property
    def theoretical(self) -> TheoreticalDistribution:
        """Return the theoretical approximation part."""
        return self._theoretical

    @property
    def pdfs(self) -> pd.DataFrame:
        """Return the theoretical approximation part."""
        if self._pdfs is None:
            self._pdfs = self._get_pdfs()
        return self._pdfs

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
