#!/usr/bin/env python3
"""Nature data set."""

from dataclasses import dataclass
from logging import debug

import pandas as pd

from data_set.data_set import DataSet


class NatureDataSet(DataSet):
    """This class represents the "Nature Data Set"."""

    @dataclass
    class Parameters(DataSet.Parameters):
        """Represent the necessary properties to load a dataset."""

        date_interval: tuple[pd.Timestamp, pd.Timestamp]

    def __init__(self, data_set_properties: Parameters) -> None:
        """Create data set without loading data."""
        super().__init__(data_set_properties)
        self._authors: pd.DataFrame = pd.DataFrame()
        self._documents: pd.DataFrame = pd.DataFrame()

    def _read_data(self) -> None:
        """Load data from the disk for further processing."""
        # the dataset contains large integers as node ids
        # we first create a table that maps these integers as strings to smaller internal node ids
        authors = pd.read_csv(self._data_set_properties.location / 'authors.csv', header=0)
        authors['internal_node_id'] = authors.index
        authors = authors.set_index('author_id')
        authors.index.name = 'external_node_id'
        self.authors = authors

        all_documents = pd.read_csv(
            self._data_set_properties.location / 'documents.csv',
            header=0,
            low_memory=False,
            index_col=0,
            parse_dates=['cover_date']
        )
        all_documents.index.name = 'document_id'
        all_documents['author_ids'] = all_documents['author_ids'].apply(eval)

        # filtering
        relevant_documents = all_documents[
            (all_documents['cover_date'] >= self._data_set_properties.date_range[0]) &
            (all_documents['cover_date'] < self._data_set_properties.date_range[1]) &
            (all_documents['subtype'].isin(['ar'])) &
            (all_documents['author_ids'].map(len) <= self._data_set_properties.max_simplex_dimension)
        ]

        def single_node_id_translator(external_node_id: int) -> int | None:
            if external_node_id in authors.index:
                return authors.loc[external_node_id, 'internal_node_id']
            return None

        def multiple_node_ids_translator(external_node_ids: list[int | None]) -> list[int]:
            translated_node_ids = list(map(single_node_id_translator, external_node_ids))
            internal_node_ids = [x for x in translated_node_ids if x is not None]
            return internal_node_ids

        relevant_documents['author_ids'] = \
            relevant_documents['author_ids'].apply(multiple_node_ids_translator)

        self.documents = relevant_documents

    def _build_simplicial_complex(self) -> None:
        """Build a simplicial complex based on the loaded data."""
        assert not self._documents.empty, 'Data is not loaded.'
        self.add_simplices(list(self.documents['author_ids']))
        debug('Simplicial complex built.')

    @property
    def authors(self) -> pd.DataFrame:
        """Get the authors of the data set."""
        return self._authors

    @authors.setter
    def authors(self, value: pd.DataFrame) -> None:
        self._authors = value

    @authors.deleter
    def authors(self) -> None:
        del self._authors

    @property
    def documents(self) -> pd.DataFrame:
        """Get the documents of the data set."""
        return self._documents

    @documents.setter
    def documents(self, value: pd.DataFrame) -> None:
        self._documents = value

    @documents.deleter
    def documents(self) -> None:
        del self._documents

    def __str__(self) -> str:
        """Return a string representation based on the data set properties."""
        date_interval = \
            f"[{self._data_set_properties.date_interval[0].strftime('%Y-%m-%d')}, " + \
            f"{self._data_set_properties.date_interval[1].strftime('%Y-%m-%d')}]"

        component_index = self._data_set_properties.component_index_from_largest
        component = f'component_{component_index}' if component_index != -1 else 'whole'
        return '\n'.join([
            'Name: Nature',
            f'Date Interval: {date_interval},'
            f'Max Dimension: {self._data_set_properties.max_dimension}',
            f'Component: {component}'
        ])
