#!/usr/bin/env python3
"""Arxiv data set."""

from dataclasses import dataclass
from logging import debug
from typing import Any

import pandas as pd

from dataset.dataset import Dataset
from dataset.arxiv_categories import ArxivField, ArxivSubCategory
from tools.logging_helper import log_function_name


class ArxivDataset(Dataset):
    """This class represents the "Nature Data Set"."""

    @dataclass
    class Parameters(Dataset.Parameters):
        """Represent the necessary properties to load a dataset."""

        date_interval: tuple[pd.Timestamp, pd.Timestamp]
        fields: list[ArxivField]
        primary_categories: list[ArxivSubCategory]

    def __init__(self, data_set_properties: Parameters) -> None:
        """Create data set without loading data."""
        super().__init__(data_set_properties)
        self._authors: pd.DataFrame = pd.DataFrame()
        self._documents: pd.DataFrame = pd.DataFrame()

    def info(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        date_interval = [
            self._data_set_properties.date_interval[0].strftime('%Y-%m-%d'),
            self._data_set_properties.date_interval[1].strftime('%Y-%m-%d'),
        ]

        if ArxivField.INVALID in self._data_set_properties.fields:
            field_names = 'Not filtered'
        else:
            field_names = [field.name for field in self._data_set_properties.fields]

        if ArxivSubCategory.INVALID in self._data_set_properties.primary_categories:
            categories = 'Not filtered'
        else:
            categories = [category.name for category in self._data_set_properties.primary_categories]

        result: dict[str, Any] = super().info()
        result.update({
            'name': 'arxiv',
            'fields': field_names,
            'categories': categories,
            'from_date': date_interval[0],
            'to_date': date_interval[1],
            'component': self._data_set_properties.component_index_from_largest,
        })
        return result

    @log_function_name
    def _read_data(self) -> None:
        """Load data from the disk for further processing."""
        # the dataset contains large integers as node ids
        # we first create a table that maps these integers as strings to smaller internal node ids
        data_file_name = self._data_set_properties.location / 'documents.csv'
        debug(f'Reading data file {str(data_file_name)}...')
        all_documents = pd.read_csv(
            data_file_name,
            index_col=0,
            parse_dates=['update_time', 'publish_time'],
            low_memory=False
        )
        debug('Reading data file done')
        all_documents['categories'] = all_documents['categories'].apply(eval)
        all_documents['authors'] = all_documents['authors'].apply(eval).apply(set).apply(list)

        # filtering
        if ArxivField.INVALID not in self._data_set_properties.fields:
            filtered_documents = all_documents[
                all_documents['field'].isin([field.value for field in self._data_set_properties.fields])
            ]

        if ArxivSubCategory.INVALID not in self._data_set_properties.primary_categories:
            filtered_documents = filtered_documents[
                filtered_documents['primary_category'].isin([
                    category.value for category in self._data_set_properties.primary_categories
                ])
            ]

        print(f'Filtered documents: {len(filtered_documents)} / ', end='')

        filtered_documents = filtered_documents[
            (filtered_documents['publish_time'] >=
                pd.Timestamp(self._data_set_properties.date_interval[0], tz='UTC')) &
            (filtered_documents['publish_time'] <
                pd.Timestamp(self._data_set_properties.date_interval[1], tz='UTC')) &
            (filtered_documents['authors'].map(len) <=
                self._data_set_properties.max_simplex_dimension)
        ]

        assert len(filtered_documents) > 0, 'No data is available with the specified properties.'

        self._documents = filtered_documents
        # self._fill_authors_table()

        print(len(filtered_documents))

    @log_function_name
    def _get_vertices_from_data(self) -> list[int]:
        """Return vertices of the loaded data."""
        return sorted(list(set(self.documents['authors'].sum())))

    @log_function_name
    def _get_interactions_from_data(self) -> list[list[int]]:
        """Build a simplicial complex based on the loaded data."""
        assert not self._documents.empty, 'Data is not loaded.'

        interactions = [
            list(authors)
            for authors in self.documents['authors']
            if len(authors) != 0
        ]

        return interactions

    @log_function_name
    def _fill_authors_table(self) -> None:

        list_of_author_ids = sorted(list(set(self.documents['authors'].sum())))
        self._authors = pd.DataFrame(list_of_author_ids, columns=['author_id'])
        self._authors.set_index('author_id', inplace=True)

        def single_author_name_translator(name: int) -> int | None:
            if name in self._authors.index:
                return self._authors.loc[name, 'author_id']
            return None

        def multiple_author_names_translator(names: list[str | None]) -> list[int]:
            author_ids = list(map(single_author_name_translator, names))
            author_ids = [x for x in author_ids if x is not None]
            return author_ids

        self.documents['authors'] = self.documents['authors'].apply(multiple_author_names_translator)

    @property
    def authors(self) -> pd.DataFrame:
        """Get the authors of the data set."""
        return self._authors

    @property
    def documents(self) -> pd.DataFrame:
        """Get the documents of the data set."""
        return self._documents
