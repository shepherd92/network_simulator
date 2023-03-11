#!/usr/bin/env python3
"""This module represents the Watts-Strogatz network model."""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

from data_set.data_set import DataSet
from model.model import Model
from network.network import Network
from network.property import BaseNetworkProperty


class WattsStrogatzModel(Model):
    """Class representing a Watts-Strogatz model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an WattsStrogatzModel."""

        edges_of_new_node: int = 0
        rewiring_probability: float = 0.

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = WattsStrogatzModel.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.Type.NUM_OF_NODES
        ))
        num_of_edges: int = data_set.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.Type.NUM_OF_EDGES
        ))
        edges_of_new_node_guess = round(num_of_edges / num_of_nodes)

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.num_nodes = num_of_nodes
        # pylint: enable=attribute-defined-outside-init
        self._parameters.edges_of_new_node = edges_of_new_node_guess

    def generate_network(self, seed: int | None = None) -> Network:
        """Build a network of the model."""
        assert isinstance(self.parameters, WattsStrogatzModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        graph: nx.Graph = nx.watts_strogatz_graph(
            self.parameters.num_nodes,
            self.parameters.edges_of_new_node,
            self.parameters.rewiring_probability,
            seed=seed
        )

        network = Network(self.parameters.max_dimension)
        network.graph = graph
        network.digraph = graph.to_directed()
        network.generate_simplicial_complex_from_graph()
        network.interactions = graph.edges

        return network

    @property
    def parameters(self) -> WattsStrogatzModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: WattsStrogatzModel.Parameters) -> None:
        self._parameters = value
