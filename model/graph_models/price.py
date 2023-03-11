#!/usr/bin/env python3
"""This module represents the Erdos-Renyi network model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from data_set.data_set import DataSet
from model.model import Model
from network.network import Network
from network.property import BaseNetworkProperty


class PriceModel(Model):
    """Class representing an Erdos-Renyi model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an PriceModel."""

        probability_degree_constant: float = 1.

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = PriceModel.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.Type.NUM_OF_NODES
        ))

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.num_nodes = num_of_nodes
        # pylint: enable=attribute-defined-outside-init

    def generate_network(self, _: int | None = None) -> Network:
        """Build a network of the model."""
        assert isinstance(self.parameters, PriceModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        node_ids = np.array(range(self.parameters.num_nodes))
        out_degrees = self._get_outdegrees()
        corrected_in_degrees = \
            np.zeros_like(out_degrees) + self.parameters.probability_degree_constant
        edges: list[list[int]] = []

        for node_id in node_ids:
            probabilities = corrected_in_degrees[:node_id] / np.sum(corrected_in_degrees[:node_id])
            if node_id == 0:
                continue
            chosen_node_ids = np.random.choice(
                range(node_id),
                replace=False,
                p=probabilities,
                size=out_degrees[node_id]
            )
            corrected_in_degrees[chosen_node_ids] += 1
            edges.extend(list(zip(list(chosen_node_ids), [node_id] * len(chosen_node_ids))))

        network = Network(self.parameters.max_dimension)
        network.graph.add_nodes_from(node_ids)
        network.graph.add_edges_from(edges)
        network.digraph = network.graph.to_directed()
        network.generate_simplicial_complex_from_graph()

        return network

    def _get_outdegrees(self):
        """Generate out degrees of each node based on a given distribution."""
        out_degrees = np.random.randint(0, 5, size=self.parameters.num_nodes)

        # The outdegree of the i-th node can be at most i.
        effective_outdegrees = \
            np.min(np.c_[out_degrees, np.array(range(self.parameters.num_nodes))], axis=1)

        return effective_outdegrees

    @property
    def parameters(self) -> PriceModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: PriceModel.Parameters) -> None:
        self._parameters = value
