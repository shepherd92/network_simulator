#!/usr/bin/env python3
"""Bianconi data set."""

from dataclasses import dataclass

import pandas as pd

from data_set.data_set import DataSet


class BianconiDataSet(DataSet):
    """This class represents the "Bianconi Data Set"."""

    @dataclass
    class Parameters(DataSet.Parameters):
        """Represent the necessary properties to load a dataset."""

    def __init__(self, data_set_properties: Parameters) -> None:
        """Create data set without loading data."""
        super().__init__(data_set_properties)
        self._data: dict[str, pd.DataFrame] = {}

    def _read_data(self) -> None:
        """Load data from the disk for further processing."""
        # the dataset contains large integers as node ids
        # we first create a table that maps these integers as strings to smaller internal node ids
        data_location = self._data_set_properties.location
        node_ids = pd.read_csv(data_location / 'nodelist_ConnComp.txt', header=None)
        node_ids['internal_node_id'] = node_ids.index
        node_ids = node_ids.set_index(0)
        node_ids.index.name = 'external_node_id'

        edges = pd.read_csv(
            data_location / 'edgelist_weighted.txt', header=None,
            names=['node_0', 'node_1', 'weight']
        )
        edges[['node_0', 'node_1']] = edges[['node_0', 'node_1']] \
            .replace(node_ids.index, node_ids['internal_node_id'])

        triangles = pd.read_csv(
            data_location / 'trilist_weighted.txt', header=None,
            names=['node_0', 'node_1', 'node_2', 'weight']
        )
        triangles[['node_0', 'node_1', 'node_2']] = \
            triangles[['node_0', 'node_1', 'node_2']] \
            .replace(node_ids.index, node_ids['internal_node_id'])

        self._data = {
            'nodes': node_ids,
            'edges': edges,
            'triangles': triangles,
        }

    def _build_simplicial_complex(self) -> None:
        """Build a simplicial complex based on the loaded data."""
        assert self._data, 'Data is not loaded.'

        for node in self._data['nodes'].values:
            self.add_simplex(node)

        for edge in self._data['edges'].values:
            self.add_simplex(edge[:2], filtration=edge[-1])

        for triangle in self._data['triangles'].values:
            self.add_simplex(triangle[:3], filtration=triangle[-1])

        # pylint: disable-next=attribute-defined-outside-init
        self._interactions = \
            self._data['nodes'].values.tolist() + \
            self._data['edges'].iloc[:, :2].values.tolist() + \
            self._data['triangles'].iloc[:, :3].values.tolist()
        pass

    def _build_graph(self) -> None:
        """Build a simple networkx graph."""
        self.generate_graph_from_simplicial_complex()

    def __str__(self) -> str:
        """Return a string representation based on the data set properties."""
        return '\n'.join([
            'Name: Bianconi',
            f'Max Dimension: {self._data_set_properties.max_dimension}',
        ])
