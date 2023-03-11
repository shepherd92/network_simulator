#!/usr/bin/env python3
"""Configuration of the Test data set."""

from pathlib import Path

from data_set.test_data_set import TestDataSet


TEST_DATA_SET_PARAMETERS = TestDataSet.Parameters(
    location=Path(),
    max_dimension=2,
    max_simplex_dimension=20,
    component_index_from_largest=-1,
)
