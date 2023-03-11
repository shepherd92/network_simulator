#!/usr/bin/env python3
"""Configuration of the Arxiv data set."""

from pathlib import Path

import pandas as pd

from data_set.arxiv_data_set import ArxivDataSet
from data_set.arxiv_categories import ArxivField, ArxivSubCategory


ARXIV_DATA_SET_PARAMETERS = ArxivDataSet.Parameters(
    location=Path('../../data/arxiv'),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    date_interval=(pd.Timestamp('2021-01-01'), pd.Timestamp('2023-12-31')),
    field=ArxivField.INVALID,
    primary_category=ArxivSubCategory.stat_ML,
)
