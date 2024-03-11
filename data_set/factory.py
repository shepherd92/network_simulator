#!/usr/bin/env python3
"""Data Loader."""

from logging import info

from config_files.data_set_config import ARXIV_DATA_SET_PARAMETERS
from config_files.data_set_config import BIANCONI_DATA_SET_PARAMETERS
from config_files.data_set_config import NATURE_DATA_SET_PARAMETERS
from config_files.data_set_config import TEST_DATA_SET_PARAMETERS

from data_set.data_set import DataSet
from data_set.arxiv_data_set import ArxivDataSet
from data_set.bianconi_data_set import BianconiDataSet
from data_set.nature_data_set import NatureDataSet
from data_set.test_data_set import TestDataSet


def load_data(data_set_to_load: DataSet.Type) -> DataSet:
    """Load the specified data set."""
    if data_set_to_load == DataSet.Type.BIANCONI:
        data = BianconiDataSet(BIANCONI_DATA_SET_PARAMETERS)
    elif data_set_to_load == DataSet.Type.NATURE:
        data = NatureDataSet(NATURE_DATA_SET_PARAMETERS)
    elif data_set_to_load == DataSet.Type.ARXIV:
        data = ArxivDataSet(ARXIV_DATA_SET_PARAMETERS)
    elif data_set_to_load == DataSet.Type.TEST:
        data = TestDataSet(TEST_DATA_SET_PARAMETERS)
    else:
        raise NotImplementedError(f'The requested data set {data_set_to_load.name} is not yet implemented')

    info(f'Data set {data_set_to_load.name} loaded.')
    return data
