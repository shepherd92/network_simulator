#!/usr/bin/env python3
"""Configuration of the Nature data set."""

from pathlib import Path

import pandas as pd

from data_set.nature_data_set import NatureDataSet


NATURE_DATA_SET_PARAMETERS = NatureDataSet.Parameters(
    location=Path('../../data/nature'),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
    date_interval=(pd.Timestamp('1900-01-01'), pd.Timestamp('2023-12-31')),
)
