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
    properties_to_calculate: list[BaseNetworkProperty.Type] = [
        BaseNetworkProperty.Type.NUM_OF_VERTICES,
        BaseNetworkProperty.Type.NUM_OF_EDGES,
        BaseNetworkProperty.Type.NUM_OF_TRIANGLES,
        BaseNetworkProperty.Type.NUM_OF_INTERACTIONS,
        # BaseNetworkProperty.Type.EDGES,
        BaseNetworkProperty.Type.AVERAGE_DEGREE,
        BaseNetworkProperty.Type.MAX_DEGREE,
        # BaseNetworkProperty.Type.AVG_CLUSTERING,
        BaseNetworkProperty.Type.NUM_OF_CONNECTED_COMPONENTS,
        BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION,
        # BaseNetworkProperty.Type.SIMPLEX_DIMENSION_DISTRIBUTION,
        # BaseNetworkProperty.Type.FACET_DIMENSION_DISTRIBUTION,
        BaseNetworkProperty.Type.DEGREE_DISTRIBUTION,
        BaseNetworkProperty.Type.VERTEX_INTERACTION_DEGREE_DISTRIBUTION,
        BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1,
        # BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2,
        BaseNetworkProperty.Type.BETTI_NUMBERS,
        # BaseNetworkProperty.Type.BETTI_NUMBERS_BY_COMPONENT,
        BaseNetworkProperty.Type.VERTICES_BY_COMPONENT,
        # BaseNetworkProperty.Type.PERSISTENCE_PAIRS,
    ]


ARXIV_DATA_SET_PARAMETERS = ArxivDataSet.Parameters(
    location=Path('../../collaboration_network/data/arxiv'),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    date_interval=(pd.Timestamp('1900-01-01'), pd.Timestamp('2024-12-31')),
    fields=[ArxivField.computer_science],
    primary_categories=[ArxivSubCategory.INVALID],  # stat_TH is an alias for math_ST
)


BIANCONI_DATA_SET_PARAMETERS = BianconiDataSet.Parameters(
    location=Path('../../data/bianconi'),
    max_dimension=3,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
)


NATURE_DATA_SET_PARAMETERS = NatureDataSet.Parameters(
    location=Path('../../data/nature'),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    date_interval=(pd.Timestamp('1900-01-01'), pd.Timestamp('2023-12-31')),
)


TEST_DATA_SET_PARAMETERS = TestDataSet.Parameters(
    location=Path(),
    max_dimension=3,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
)
