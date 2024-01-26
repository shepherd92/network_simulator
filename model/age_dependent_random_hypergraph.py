#!/usr/bin/env python3
"""This module represents the age-dependent random hypergraph network model."""

from __future__ import annotations

from dataclasses import dataclass
from logging import info

import numpy as np
import numpy.typing as npt

from data_set.data_set import DataSet
from model.model import Model
from network.finite_network import FiniteNetwork
from network.property import BaseNetworkProperty


class AgeDependentRandomHypergraphModel(Model):
    """Class representing an age-dependent random hypergraph network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an AgeDependentRandomHypergraphModel."""

        num_of_nodes: int = 0  # actual number of nodes
        num_of_interactions: int = 0  # number of interactions in the hypergraph
        beta: float = 0.1
        gamma: float = 0.5
        gamma_prime: float = 0.5
        torus_dimension: int = 0
        torus_size_in_1_dimension: float = 0.5  # the total size of the torus

        def to_numpy(self) -> npt.NDArray[np.float_]:
            """Return the parameters as a numpy array."""
            return np.array([
                self.max_dimension,
                self.network_size,
                self.num_of_nodes,
                self.num_of_interactions,
                self.beta,
                self.gamma,
                self.gamma_prime,
                self.torus_dimension,
                self.torus_size_in_1_dimension,
            ])

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        super().__init__()
        self._parameters = AgeDependentRandomHypergraphModel.Parameters()
        self._random_number_generator = np.random.default_rng()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.network_size = num_of_nodes
        # pylint: enable=attribute-defined-outside-init

        print('\n'.join([
            '\nADRCM model paramerers after setting from data set:',
            f'size        = {self._parameters.network_size}',
            f'max_dim     = {self._parameters.max_dimension}',
            f'beta        = {self._parameters.beta:4f}',
            f'gamma       = {self._parameters.gamma:4f}',
            f'gamma_prime = {self._parameters.gamma_prime:4f}\n',
        ]))

    def generate_finite_network(self, seed: int) -> FiniteNetwork:
        """Build a network of the model."""
        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed}.')
        assert isinstance(self.parameters, AgeDependentRandomHypergraphModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        self._random_number_generator = np.random.default_rng(seed)
        self.parameters.torus_size_in_1_dimension = \
            self.parameters.network_size ** (1. / self.parameters.torus_dimension)

        self.parameters.num_of_nodes = self._random_number_generator.poisson(lam=self.parameters.network_size)
        self.parameters.num_of_interactions = self._random_number_generator.poisson(lam=self.parameters.network_size)

        node_birth_times = self._random_number_generator.random(size=self.parameters.num_of_nodes)
        node_positions = self._random_number_generator.uniform(
            -self.parameters.torus_size_in_1_dimension / 2,
            +self.parameters.torus_size_in_1_dimension / 2,
            size=(self.parameters.num_of_nodes, self.parameters.torus_dimension)
        )

        interaction_birth_times = self._random_number_generator.random(size=self.parameters.num_of_interactions)
        interaction_positions = self._random_number_generator.uniform(
            -self.parameters.torus_size_in_1_dimension / 2,
            +self.parameters.torus_size_in_1_dimension / 2,
            size=(self.parameters.num_of_interactions, self.parameters.torus_dimension)
        )

        network = FiniteNetwork(self.parameters.max_dimension)
        network._interactions = self.generate_finite_network_interactions(
            node_birth_times, node_positions, interaction_birth_times, interaction_positions
        )

        node_ids = np.array(range(self.parameters.num_of_nodes))[:, np.newaxis]
        interaction_ids = np.array(range(self.parameters.num_of_interactions))[:, np.newaxis]
        if self.parameters.torus_dimension == 1:
            network.vertex_positions = {
                node_id: (node_position, node_birth_time)
                for node_id, node_position, node_birth_time
                in np.c_[node_ids, node_positions, node_birth_times]
            }
            network.interaction_positions = {
                interaction_id: (interaction_position, interaction_birth_time)
                for interaction_id, interaction_position, interaction_birth_time
                in np.c_[interaction_ids, interaction_positions, interaction_birth_times]
            }
        elif self.parameters.torus_dimension == 2:
            network.vertex_positions = {
                node_id: (node_position_0, node_position_1, node_birth_time)
                for node_id, node_position_0, node_position_1, node_birth_time
                in np.c_[node_ids, node_positions, node_birth_times]
            }
            network.interaction_positions = {
                interaction_id: (interaction_position_0, interaction_position_1, interaction_birth_time)
                for interaction_id, interaction_position_0, interaction_position_1, interaction_birth_time
                in np.c_[interaction_ids, interaction_positions, interaction_birth_times]
            }

        network.add_simplices_batch(node_ids)
        network.add_simplices(network._interactions)

        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed} done.')
        return network

    def generate_finite_network_interactions(
        self,
        node_birth_times: npt.NDArray[np.float_],
        node_positions: npt.NDArray[np.float_],
        interaction_birth_times: npt.NDArray[np.float_],
        interaction_positions: npt.NDArray[np.float_],
    ) -> list[list[int]]:
        """Generate the connections of a finite network."""

        # interactions: dict[int, list[npt.NDArray[np.int_]]] = {}
        interactions: list[list[int]] = []
        for interaction_birth_time, interaction_position in zip(
            interaction_birth_times.tolist(),
            interaction_positions.tolist(),
        ):
            distances = self._distances(
                np.c_[np.full(self.parameters.num_of_nodes, interaction_position), node_positions]
            )
            distances_d = distances**self.parameters.torus_dimension
            profile_function_argument = distances_d / self.parameters.beta \
                * node_birth_times**(self.parameters.gamma) \
                * interaction_birth_time**(self.parameters.gamma_prime)
            interaction_members = np.where(profile_function_argument < 0.5)[0]
            interactions.append(interaction_members.tolist())

        return interactions

    def _distances(self, position_pairs: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:

        dimensions = position_pairs.shape[1] // 2

        # calculate the distances by dimension inside the boundaries without going around the edges
        distances_by_dimension_inside = np.abs(position_pairs[:, :dimensions] - position_pairs[:, dimensions:])

        # if the distance in one dimension is greater than the half of the size of the torus,
        # then it is better to go around the edge, when the distance in that direction is
        # 1 - distance_inside
        distances_by_dimension = np.where(
            distances_by_dimension_inside < 0.5 * self.parameters.torus_size_in_1_dimension,
            distances_by_dimension_inside,
            self.parameters.torus_size_in_1_dimension - distances_by_dimension_inside
        )

        # the final distances are the length of the difference vectors
        distances = np.linalg.norm(distances_by_dimension, axis=1)

        return distances

    def set_model_parameters_from_tuple(self, parameters_tuple: tuple[int]) -> None:
        """Convert a tuple to ModelParamters. Used for model optimization."""
        self.parameters = AgeDependentRandomHypergraphModel.Parameters(*parameters_tuple)

    def get_info_as_dict(self) -> dict[str, int | float]:
        """Return a dict representation based on the model properties."""
        return {
            'max_dimension': self.parameters.max_dimension,
            'network_size': self.parameters.network_size,
            'num_of_nodes': self.parameters.num_of_nodes,
            'parameter_beta': self.parameters.beta,
            'parameter_gamma': self.parameters.gamma,
            'parameter_gamma_prime': self.parameters.gamma_prime,
            'torus_dimension': self.parameters.torus_dimension,
            'torus_size_in_1_dimension': self.parameters.torus_size_in_1_dimension,
        }

    @property
    def parameters(self) -> AgeDependentRandomHypergraphModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: AgeDependentRandomHypergraphModel.Parameters) -> None:
        self._parameters = value
