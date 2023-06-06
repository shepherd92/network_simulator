#!/usr/bin/env python3
"""Configuration of the Arxiv data set."""

from pathlib import Path

import pandas as pd

from data_set.arxiv_data_set import ArxivDataSet
from data_set.arxiv_categories import ArxivField, ArxivSubCategory


ARXIV_DATA_SET_PARAMETERS = ArxivDataSet.Parameters(
    location=Path('../../data/arxiv'),
    max_dimension=3,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    date_interval=(pd.Timestamp('1900-01-01'), pd.Timestamp('2023-12-31')),
    fields=[ArxivField.engineering],
    primary_categories=[ArxivSubCategory.INVALID],
)

# stat_TH is an alias for math_ST
