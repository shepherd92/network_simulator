#!/usr/bin/env python3
"""Arxiv data set."""

from dataclasses import dataclass
from logging import debug

import networkx as nx
import pandas as pd

from data_set.data_set import DataSet
from data_set.arxiv_categories import ArxivField, ArxivSubCategory


class ArxivDataSet(DataSet):
    """This class represents the "Nature Data Set"."""

    @dataclass
    class Parameters(DataSet.Parameters):
        """Represent the necessary properties to load a dataset."""

        date_interval: tuple[pd.Timestamp, pd.Timestamp]
        fields: list[ArxivField]
        primary_categories: list[ArxivSubCategory]

    def __init__(self, data_set_properties: Parameters) -> None:
        """Create data set without loading data."""
        super().__init__(data_set_properties)
        self._data_set_properties = data_set_properties
        self._authors: pd.DataFrame = pd.DataFrame()
        self._documents: pd.DataFrame = pd.DataFrame()

    def _read_data(self) -> None:
        """Load data from the disk for further processing."""
        # the dataset contains large integers as node ids
        # we first create a table that maps these integers as strings to smaller internal node ids
        debug('Reading data file...')
        all_documents = pd.read_csv(
            self._data_set_properties.location / 'documents.csv',
            index_col=0,
            parse_dates=['update_time', 'publish_time'],
            low_memory=False
        )
        debug('done')
        all_documents['categories'] = all_documents['categories'].apply(eval)
        all_documents['authors'] = all_documents['authors'].apply(eval)

        # filtering
        filtered_documents = all_documents[
            (all_documents['publish_time'] >=
                pd.Timestamp(self._data_set_properties.date_interval[0], tz='UTC')) &
            (all_documents['publish_time'] <
                pd.Timestamp(self._data_set_properties.date_interval[1], tz='UTC')) &
            (all_documents['authors'].map(len) <= self._data_set_properties.max_simplex_dimension)
        ]

        if ArxivField.INVALID not in self._data_set_properties.fields:
            filtered_documents = filtered_documents[
                filtered_documents['field'].isin([field.value for field in self._data_set_properties.fields])
            ]

        if ArxivSubCategory.INVALID not in self._data_set_properties.primary_categories:
            filtered_documents = filtered_documents[
                filtered_documents['primary_category'].isin([
                    category.value for category in self._data_set_properties.primary_categories
                ])
            ]

        assert len(filtered_documents) > 0, 'No data is available with the specified properties.'

        self._documents = filtered_documents
        # self._fill_authors_table()

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

    def _build_simplicial_complex(self) -> None:
        """Build a simplicial complex based on the loaded data."""
        assert not self._documents.empty, 'Data is not loaded.'

        debug('Building simplicial complex...')

        self.add_simplices(self.documents['authors'])

        # pylint: disable-next=attribute-defined-outside-init
        self._interactions = [
            list(simplex)
            for simplex in self.documents['authors']
            if len(simplex) != 0
        ]

        debug('done.')

    def _build_graph(self) -> None:
        """Build a simple networkx graph."""
        self.generate_graph_from_simplicial_complex()

    @property
    def authors(self) -> pd.DataFrame:
        """Get the authors of the data set."""
        return self._authors

    @property
    def documents(self) -> pd.DataFrame:
        """Get the documents of the data set."""
        return self._documents

    def __str__(self) -> str:
        """Return a string representation based on the data set properties."""
        date_interval = \
            f"[{self._data_set_properties.date_interval[0].strftime('%Y-%m-%d')}, " + \
            f"{self._data_set_properties.date_interval[1].strftime('%Y-%m-%d')}]"

        component_index = self._data_set_properties.component_index_from_largest
        component = f'component_{component_index}' if component_index != -1 else 'whole'

        if ArxivField.INVALID in self._data_set_properties.fields:
            field_names = 'Not filtered'
        else:
            field_names = [field.name for field in self._data_set_properties.fields]

        if ArxivSubCategory.INVALID in self._data_set_properties.primary_categories:
            categories = 'Not filtered'
        else:
            categories = [category.name for category in self._data_set_properties.primary_categories]

        return '\n'.join([
            'Name: Arxiv',
            f'Fields: {field_names}',
            f'Category: {categories}',
            f'Number of vertices: {self._simplicial_complex.num_vertices()}',
            f'Number of documents: {len(self._interactions)}',
            f'Number of simplices: {self._simplicial_complex.num_simplices()}',
            f'Date Interval: {date_interval},',
            f'Max Dimension: {self._data_set_properties.max_dimension}',
            f'Number of components: {len(list(nx.connected_components(self._graph)))}',
            f'Component: {component}',
        ])
