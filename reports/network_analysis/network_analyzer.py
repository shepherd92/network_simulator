#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import (
    create_power_law_fitting_parameters,
)
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from reports.plotting_helper import (
    PaddingSide,
    check_calculated,
    plot_approximation_value_counts_log,
    plot_value_counts,
)

class NetworkAnalyzer:
    """Class to analyze finite hypergraph models."""

    @dataclass
    class Parameters:
        """Parameters for the network analysis."""
        save_directory: Path
        power_law_fitting_minimum_value : float

    def __init__(self, parameters: Parameters):
        """Initialize the network analyzer with given parameters."""
        self._parameters = parameters

    @check_calculated('Vertex--edge degree distribution')
    def report_vertex_edge_degree_distribution(
        self,
        empirical_distribution: EmpiricalDistribution,
        axes: plt.Axes,
        **kwargs,
    ) -> None:

        approximation = DistributionApproximation(
            empirical_distribution,
            TheoreticalDistribution.Type.POWER_LAW
        )
        fitting_parameters = create_power_law_fitting_parameters(self.power_law_fitting_minimum_value)
        approximation.fit(fitting_parameters)

        approximation.save(self.save_directory / 'total_degree_distribution')
        plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


    @check_calculated('Edge--triangle degree distribution')
    def report_edge_triangle_degree_distribution(
        self,
        empirical_distribution: EmpiricalDistribution,
        axes: plt.Axes,
    ) -> None:

        approximation = DistributionApproximation(
            empirical_distribution,
            TheoreticalDistribution.Type.POWER_LAW
        )
        fitting_parameters = create_power_law_fitting_parameters(self.power_law_fitting_minimum_value)
        approximation.fit(fitting_parameters)

        approximation.save(self.save_directory / 'ho_degree_distribution_1')
        plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


    @check_calculated('Triangle--tetrahedron degree distribution')
    def report_triangle_tetrahedron_degree_distribution(
        self,
        empirical_distribution: EmpiricalDistribution,
        axes: plt.Axes,
    ) -> None:
        """Report triangle--tetrahedron degree distribution."""
        approximation = DistributionApproximation(
            empirical_distribution,
            TheoreticalDistribution.Type.POWER_LAW
        )
        fitting_parameters = create_power_law_fitting_parameters(self.power_law_fitting_minimum_value)
        approximation.fit(fitting_parameters)

        approximation.save(self.save_directory / 'ho_degree_distribution_2')
        plot_approximation_value_counts_log(approximation, np.nan, PaddingSide.NONE, axes)


    @check_calculated('Simplex dimension distribution')
    def report_simplex_dimension_distribution(
        self,
        empirical_distribution: EmpiricalDistribution,
        axes: plt.Axes,
    ) -> None:
        """Report simplex dimension distribution."""
        approximation = DistributionApproximation(
            empirical_distribution,
            TheoreticalDistribution.Type.POWER_LAW
        )
        fitting_parameters = create_power_law_fitting_parameters(self.power_law_fitting_minimum_value)
        approximation.fit(fitting_parameters)

        approximation.save(self.save_directory / 'simplex_dimension_distribution')
        plot_value_counts(empirical_distribution.calc_value_counts(), axes)

    @property
    def save_directory(self) -> Path:
        """Get the directory where reports are saved."""
        return self._parameters.save_directory

    @property
    def power_law_fitting_minimum_value(self) -> float:
        """Return the power law fitting minimum value."""
        return self._parameters.power_law_fitting_minimum_value
