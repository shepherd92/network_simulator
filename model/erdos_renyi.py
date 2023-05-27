#!/usr/bin/env python3
"""Erdos-Renyi network model."""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

from data_set.data_set import DataSet
from model.model import Model
from network.finite_network import FiniteNetwork
from network.property import BaseNetworkProperty


class ErdosRenyiModel(Model):
    """Erdos-Renyi network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an ErdosRenyiModel."""

        edge_probability: float = 0.

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = ErdosRenyiModel.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.Type.NUM_OF_NODES
        ))
        num_of_edges: int = data_set.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.Type.NUM_OF_EDGES
        ))
        edge_probability_guess = num_of_edges / (num_of_nodes * (num_of_nodes - 1) / 2)

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.num_nodes = num_of_nodes
        # pylint: enable=attribute-defined-outside-init
        self._parameters.edge_probability = edge_probability_guess

    def generate_finite_network(self, seed: int | None = None) -> FiniteNetwork:
        """Build a network of the model."""
        graph: nx.Graph = nx.erdos_renyi_graph(
            self._parameters.num_nodes,
            self._parameters.edge_probability,
            seed=seed
        )

        network = FiniteNetwork(self._parameters.max_dimension)
        network.graph = graph
        network.digraph = graph.to_directed()
        network.generate_simplicial_complex_from_graph()
        network._interactions = graph.edges
        network._facets = graph.edges

        return network

    @property
    def parameters(self) -> ErdosRenyiModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: ErdosRenyiModel.Parameters) -> None:
        self._parameters = value
