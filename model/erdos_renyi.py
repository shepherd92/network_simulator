#!/usr/bin/env python3
"""Erdos-Renyi network model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import NewType

import networkx as nx

from dataset.dataset import Dataset
from model.model import Model
from network.finite_clique_complex import FiniteCliqueComplex
from network.property import BaseNetworkProperty


InfiniteNetworkSet = NewType('InfiniteNetworkSet', None)


class ErdosRenyiModel(Model):
    """Erdos-Renyi network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an ErdosRenyiModel."""

        edge_probability: float = 0.

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = ErdosRenyiModel.Parameters()

    def set_relevant_parameters_from_dataset(self, dataset: Dataset) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = dataset.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.num_of_vertices
        ))
        num_of_edges: int = dataset.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.num_of_edges
        ))
        edge_probability_guess = num_of_edges / (num_of_nodes * (num_of_nodes - 1) / 2)

        self._parameters = ErdosRenyiModel.Parameters(
            max_dimension=dataset.max_dimension,
            network_size=num_of_nodes,
            edge_probability=edge_probability_guess,
        )

        print(self)

    def generate_finite_network(self, seed: int | None = None) -> FiniteCliqueComplex:
        """Build a network of the model."""
        graph: nx.Graph = nx.erdos_renyi_graph(
            self._parameters.network_size,
            self._parameters.edge_probability,
            seed=seed
        )

        network = FiniteCliqueComplex(self._parameters.max_dimension)
        network.graph = graph
        network.generate_simplicial_complex_from_graph()
        network.interactions = graph.edges

        return network

    def generate_infinite_network_set(self, num_of_networks: int, seed: int) -> InfiniteNetworkSet:
        raise NotImplementedError

    @property
    def parameters(self) -> ErdosRenyiModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: ErdosRenyiModel.Parameters) -> None:
        self._parameters = value
