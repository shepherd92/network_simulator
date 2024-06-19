#!/usr/bin/env python3
"""Configuration of the data sets."""

from pathlib import Path
from typing import NamedTuple

import pandas as pd

from data_set.data_set import DataSet
from data_set.arxiv_data_set import ArxivDataSet
from data_set.arxiv_categories import ArxivField, ArxivSubCategory
from data_set.bianconi_data_set import BianconiDataSet
from data_set.nature_data_set import NatureDataSet
from data_set.test_data_set import TestDataSet
from network.property import BaseNetworkProperty


class DataSetConfig(NamedTuple):
    """Data set configuration."""

    type_: DataSet.Type = DataSet.Type.ARXIV
    plot: bool = True
    properties_to_calculate: list[BaseNetworkProperty] = [
        BaseNetworkProperty.num_of_vertices,
        BaseNetworkProperty.num_of_edges,
        BaseNetworkProperty.num_of_triangles,
        BaseNetworkProperty.num_of_interactions,
        # BaseNetworkProperty.edges,
        BaseNetworkProperty.mean_degree,
        BaseNetworkProperty.max_degree,
        BaseNetworkProperty.mean_clustering,
        BaseNetworkProperty.num_of_connected_components,
        BaseNetworkProperty.interaction_vertex_degree_distribution,
        BaseNetworkProperty.simplex_dimension_distribution,
        # BaseNetworkProperty.facet_dimension_distribution,
        # BaseNetworkProperty.in_degree_distribution,
        # BaseNetworkProperty.out_degree_distribution,
        BaseNetworkProperty.vertex_interaction_degree_distribution,
        BaseNetworkProperty.edge_interaction_degree_distribution,
        BaseNetworkProperty.vertex_edge_degree_distribution,
        BaseNetworkProperty.edge_triangle_degree_distribution,
        # BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
        BaseNetworkProperty.betti_numbers,
        # BaseNetworkProperty.betti_numbers_by_component,
        BaseNetworkProperty.num_of_vertices_by_component,
        BaseNetworkProperty.persistence_intervals,
        # BaseNetworkProperty.persistence_pairs,
    ]


ARXIV_DATA_SET_PARAMETERS = ArxivDataSet.Parameters(
    location=Path('../../collaboration_network/data/arxiv'),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    weighted=True,
    date_interval=(pd.Timestamp('1900-01-01'), pd.Timestamp('2025-01-01')),
    fields=[ArxivField.mathematics],
    primary_categories=[ArxivSubCategory.INVALID],  # stat_TH is an alias for math_ST
)


TEST_DATA_SET_PARAMETERS = TestDataSet.Parameters(
    location=Path(),
    max_dimension=3,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    weighted=True,
)


BIANCONI_DATA_SET_PARAMETERS = BianconiDataSet.Parameters(
    location=Path('../../data/bianconi'),
    max_dimension=3,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    weighted=False,
)


NATURE_DATA_SET_PARAMETERS = NatureDataSet.Parameters(
    location=Path('../../data/nature'),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    weighted=False,
    date_interval=(pd.Timestamp('1900-01-01'), pd.Timestamp('2023-12-31')),
)
