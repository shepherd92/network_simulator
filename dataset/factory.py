#!/usr/bin/env python3
"""Data Loader."""

from logging import info

from config_files.dataset_config import DatasetConfig

from dataset.dataset import Dataset
from dataset.arxiv_dataset import ArxivDataset
from dataset.test_dataset import TestDataset


def load_data(data_set_to_load: Dataset.Type) -> Dataset:
    """Load the specified data set."""
    if data_set_to_load == Dataset.Type.ARXIV:
        data = ArxivDataset(DatasetConfig.arxiv_dataset_parameters)
    elif data_set_to_load == Dataset.Type.TEST:
        data = TestDataset(DatasetConfig.test_dataset_parameters)
    else:
        raise NotImplementedError(f'The requested data set {data_set_to_load.name} is not yet implemented')

    info(f'Data set {data_set_to_load.name} loaded.')
    return data
