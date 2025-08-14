#!/usr/bin/env python3
"""Configuration of the data sets."""

from pathlib import Path
from typing import NamedTuple

import pandas as pd

from dataset.dataset import Dataset
from dataset.arxiv_dataset import ArxivDataset
from dataset.arxiv_categories import ArxivField, ArxivSubCategory
from dataset.test_dataset import TestDataset
from network.property import BaseNetworkProperty


class DatasetConfig(NamedTuple):
    """Data set configuration."""

    type_: Dataset.Type = Dataset.Type.ARXIV

    # plotting
    plot: bool = True
    plot_entire_network: bool = False
    plot_network_giant_component: bool = True
    power_law_fitting_minimum_value: float = 10.

    properties_to_calculate: list[BaseNetworkProperty] = [
        BaseNetworkProperty.num_of_vertices,
        BaseNetworkProperty.num_of_edges,
        BaseNetworkProperty.num_of_triangles,
        BaseNetworkProperty.num_of_interactions,
        BaseNetworkProperty.mean_degree,
        BaseNetworkProperty.max_degree,
        BaseNetworkProperty.interaction_vertex_degree_distribution,
        BaseNetworkProperty.simplex_dimension_distribution,
        BaseNetworkProperty.vertex_interaction_degree_distribution,
        BaseNetworkProperty.edge_interaction_degree_distribution,
        BaseNetworkProperty.vertex_edge_degree_distribution,
        BaseNetworkProperty.edge_triangle_degree_distribution,
        BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
        BaseNetworkProperty.betti_numbers,
        BaseNetworkProperty.betti_numbers_by_component,
        BaseNetworkProperty.num_of_vertices_by_component,
        # BaseNetworkProperty.persistence_intervals,
        # BaseNetworkProperty.persistence_pairs,
    ]

    arxiv_dataset_parameters = ArxivDataset.Parameters(
        location=Path('../../collaboration_network/data/arxiv'),
        max_dimension=2,
        max_simplex_dimension=20,
        component_index_from_largest=-1,
        weighted=False,
        date_interval=(pd.Timestamp('1900-01-01'), pd.Timestamp('2026-01-01')),
        fields=[ArxivField.statistics],
        primary_categories=[ArxivSubCategory.INVALID],  # stat_TH is an alias for math_ST
    )

    test_dataset_parameters = TestDataset.Parameters(
        location=Path(),
        max_dimension=3,
        max_simplex_dimension=20,
        component_index_from_largest=-1,
        weighted=True,
    )
