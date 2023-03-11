#!/usr/bin/env python3
"""Configuration of the Bianconi data set."""

from pathlib import Path

from data_set.bianconi_data_set import BianconiDataSet


BIANCONI_DATA_SET_PARAMETERS = BianconiDataSet.Parameters(
    location=Path('../../data/bianconi'),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
)
