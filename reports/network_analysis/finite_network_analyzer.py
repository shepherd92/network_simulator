#!/usr/bin/env python3
"""This module is responsible for analyzing a simplicial complex."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd

from reports.network_analysis.network_analyzer import NetworkAnalyzer
from reports.plotting_helper import (
    check_calculated,
    plot_value_counts,
    plot_value_counts_log,
)


class FiniteNetworkAnalyzer(NetworkAnalyzer):
    """Class to analyze finite hypergraph models."""

    @dataclass
    class Parameters(NetworkAnalyzer.Parameters):
        """Parameters for the network analysis."""
        plot_entire_network: bool
        plot_network_giant_component: bool
        plot_network_determined_positions: bool

    def __init__(self, parameters: Parameters):
        """Initialize the network analyzer with given parameters."""
        self._parameters = parameters

    @check_calculated('Betti numbers')
    def report_betti_numbers(
        self,
        betti_numbers: list[int],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report Betti numbers."""
        betti_numbers_array = np.c_[np.arange(len(betti_numbers)), np.array(betti_numbers)]
        plot_value_counts(betti_numbers_array, axes)
        pd.DataFrame(
            betti_numbers_array,
            columns=['dimension', 'betti_number'],
            dtype=np.int32,
        ).to_csv(self.save_directory / 'betti_numbers.csv', index=False)

    @check_calculated('Betti number 1 by component')
    def report_betti_number_1_by_component(
        self,
        betti_numbers_by_component: list[list[int]],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report Betti number 1 by component."""
        betti_numbers_1 = np.array([
            betti_numbers_in_component[1]
            for betti_numbers_in_component in betti_numbers_by_component
        ])

        values_to_plot = np.c_[np.arange(len(betti_numbers_1)), betti_numbers_1,]
        plot_value_counts_log(values_to_plot, axes)
        data_frame = pd.DataFrame(
            betti_numbers_by_component,
            columns=list(range(len(betti_numbers_by_component[0]))),
            dtype=np.int32,
        )
        data_frame.index.name = 'component_index'
        data_frame.to_csv(self.save_directory / 'betti_numbers_by_component.csv')

    @check_calculated('Vertices by component')
    def report_vertices_by_component(
        self,
        vertices_by_component: npt.NDArray[np.int_],
        axes: plt.Axes,
        **_,
    ) -> None:
        """Report vertices by component."""
        values_to_plot = np.c_[
            np.arange(1, len(vertices_by_component) + 1),
            vertices_by_component,
        ]
        plot_value_counts_log(values_to_plot, axes)
        pd.DataFrame(
            values_to_plot,
            columns=['component_index', 'num_of_vertices'],
            dtype=np.int32,
        ).to_csv(self.save_directory / 'vertices_by_component.csv', index=False)
