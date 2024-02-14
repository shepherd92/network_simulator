#!/usr/bin/env python3
"""Theoretical approximation of an empirical distribution."""

from pathlib import Path
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd

from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import (
    create_fitting_parameters_power_law_data_set,
    create_theoretical_distribution
)
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class DistributionApproximation:
    """Encloses an empirical degree distribution and its approximation."""

    class TestResult(NamedTuple):
        """Results of the statistical tests for a scalar network property."""

        probability_plot_r_value: float = np.nan
        qq_plot_r_value: float = np.nan
        point_value: float = np.nan
        point_p_value: float = np.nan
        kolmogorov_smirnov: float = np.nan

        def save(self, save_dir: Path) -> None:
            """Save the main parameters to the given file as a pandas data frame."""
            info = self.get_info_as_dict()
            data_frame = pd.DataFrame(info, index=[0])
            data_frame.to_csv(save_dir / 'test_results.csv', index=False)

        def get_info_as_dict(self) -> dict[str, int | float]:
            """Return a dict representation based on the distribution properties."""
            info_dict: dict[str, int | float] = {
                'probability_plot_r_value': self.probability_plot_r_value,
                'qq_plot_r_value': self.qq_plot_r_value,
                'kolmogorov_smirnov_statistic': self.kolmogorov_smirnov,
                'point_value': self.point_value,
                'point_p_value': self.point_p_value,
            }

            return info_dict

        def __str__(self) -> str:
            """Return string representation for reporting."""
            result = '\n'.join([
                f'probability_plot_r_value: {self.probability_plot_r_value:.4f}',
                f'qq_plot_r_value: {self.qq_plot_r_value:.4f}',
                f'Kolmogorov-Smirnov statistic: {self.kolmogorov_smirnov:.4f}',
                f'point value: {self.point_value:.4f}',
                f'point p-value: {self.point_p_value:.4f}',
            ])
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
            kolmogorov_smirnov=kolmogorov_smirnov,
            point_value=point_value,
            point_p_value=p_value,
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

    def generate_qq_plot_points(self, normalize: bool = True) -> npt.NDArray[np.float_]:
        """Return the points of the probability points."""
        number_of_values = len(self.empirical.value_sequence)
        if number_of_values == 0:
            return np.empty([0, 2])
        quantiles_to_calculate = np.linspace(1. / number_of_values, 1., number_of_values, endpoint=True)
        theoretical_quantiles = self.theoretical.calc_quantiles(quantiles_to_calculate)
        empirical_quantiles = self.empirical.value_sequence
        points = np.c_[theoretical_quantiles, empirical_quantiles]
        finite_points = points[np.isfinite(points).all(axis=1)]
        mean = self.empirical.value_sequence.mean()
        std = self.empirical.value_sequence.std()
        return (finite_points - mean) / std if normalize else finite_points

    def save(self, save_directory: Path) -> None:
        """Save all information of the approximation."""
        save_directory.mkdir(parents=True, exist_ok=True)
        self.save_info(save_directory / 'distribution_info.csv')

        if len(self.empirical.value_sequence) == 0:
            return

        self.pdfs.to_csv(save_directory / 'pdfs.csv', float_format='%.9f')
        np.savetxt(save_directory / 'value_sequence.csv', self.empirical.value_sequence, delimiter=',')
        if self.empirical.value_sequence.dtype == np.int_:
            value_counts = self.empirical.calc_value_counts()
            np.savetxt(save_directory / 'value_counts.csv', value_counts, delimiter=',')

        self.empirical.save_histogram(
            EmpiricalDistribution.HistogramType.LINEAR,
            save_directory / 'histogram_linear.csv'
        )
        self.empirical.save_histogram(
            EmpiricalDistribution.HistogramType.LOGARITHMIC,
            save_directory / 'histogram_logarithmic.csv'
        )

        confidence_levels = [0.90, 0.95, 0.99]
        confidence_intervals = self.get_confidence_intervals(confidence_levels)
        confidence_intervals.to_csv(save_directory / 'confidence_intervals.csv', float_format='%.4f')

        quantiles_to_calculate = [0.00, 0.25, 0.50, 0.75, 1.00]
        quantiles = self.get_quantiles(quantiles_to_calculate)
        quantiles.to_csv(save_directory / 'quantiles.csv', float_format='%.4f')

        probability_plot_points = self.generate_probability_plot_points()
        pd.DataFrame(
            probability_plot_points,
            columns=['theoretical', 'empirical']
        ).to_csv(save_directory / 'probability_plot_points.csv', float_format='%.4f', index=False)

        qq_plot_points = self.generate_qq_plot_points()
        pd.DataFrame(
            qq_plot_points,
            columns=['theoretical', 'empirical']
        ).to_csv(save_directory / 'qq_plot_points.csv', float_format='%.4f', index=False)

    def _get_pdfs(self) -> pd.DataFrame:
        """Get the pdfs of the two distributions."""
        num_of_points = 1000
        x_values = np.linspace(
            self.empirical.domain.min_ - 1.0 * self.empirical.domain.length,
            self.empirical.domain.max_ + 1.0 * self.empirical.domain.length,
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

    def get_quantiles(self, quantiles_to_calculate: list[float]) -> pd.DataFrame:
        """Get confidence intervals for the specified confidence levels."""
        theoretical_quantiles = self.theoretical.calc_quantiles(np.array(quantiles_to_calculate))
        empirical_quantiles = self.empirical.calc_quantiles(np.array(quantiles_to_calculate))

        index = [int(quantile*100) for quantile in quantiles_to_calculate]
        quantiles = pd.DataFrame(
            np.c_[empirical_quantiles, theoretical_quantiles],
            columns=['empirical', 'theoretical'],
            index=index
        )
        return quantiles

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

        return joint_info

    @staticmethod
    def _calc_r_squared_value(points: npt.NDArray) -> float:

        ss_res = np.sum((points[:, 1] - points[:, 0])**2)
        ss_tot = np.sum((points[:, 1] - points[:, 1].mean())**2)
        r_value = 1 - ss_res / ss_tot
        return r_value

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


def guess_power_law_exponent(distribution: EmpiricalDistribution) -> float:
    """Guess the value of parameter gamma."""
    approximation = DistributionApproximation(
        distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_parameters = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_parameters)
    return approximation.theoretical.parameters.exponent
