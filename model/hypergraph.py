#!/usr/bin/env python3
"""This module represents the age-dependent random hypergraph network model."""

from __future__ import annotations

from dataclasses import dataclass
from logging import info

import numpy as np
import numpy.typing as npt
import pandas as pd
from tqdm import tqdm

# pylint: disable-next=no-name-in-module
from cpp_modules.build.hypergraph import (  # type: ignore
    generate_finite_network_cpp,
    generate_infinite_networks_cpp,
)
from data_set.data_set import DataSet
from distribution.approximation import guess_power_law_exponent
from distribution.empirical_distribution import EmpiricalDistribution
from model.model import Model
from network.finite_network import FiniteNetwork
from network.infinite_network import InfiniteNetwork, InfiniteNetworkSet
from network.property import BaseNetworkProperty


class HypergraphModel(Model):
    """Class representing an age-dependent random hypergraph network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an AgeDependentRandomHypergraphModel."""

        interaction_intensity: float = 0.  # expected number of interactions = network_size * interaction_intensity
        beta: float = 0.1
        gamma: float = 0.5
        gamma_prime: float = 0.5
        torus_dimension: int = 0
        torus_size_in_1_dimension: float = 0.5  # the total size of the torus

        def to_numpy(self) -> npt.NDArray[np.float_]:
            """Return the parameters as a numpy array."""
            return np.array([
                self.network_size,
                self.interaction_intensity,
                self.beta,
                self.gamma,
                self.gamma_prime,
            ])

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = HypergraphModel.Parameters()
        self._random_number_generator = np.random.default_rng()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_VERTICES)
        num_of_interactions: int = data_set.calc_base_property(
            BaseNetworkProperty.Type.NUM_OF_INTERACTIONS
        )
        interaction_intensity = num_of_interactions / num_of_nodes

        vertex_interaction_degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.Type.VERTEX_INTERACTION_DEGREE_DISTRIBUTION
        )
        vertex_interaction_exponent = guess_power_law_exponent(vertex_interaction_degree_distribution)
        gamma_guess = 1. / (vertex_interaction_exponent - 1.)
        assert gamma_guess > 0. and gamma_guess < 1., f'gamma_guess = {gamma_guess}'
        document_degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION
        )
        interaction_vertex_exponent = guess_power_law_exponent(document_degree_distribution)
        gamma_prime_guess = 1. / (interaction_vertex_exponent - 1.)
        assert gamma_prime_guess > 0. and gamma_prime_guess < 1., f'gamma_prime_guess = {gamma_prime_guess}'

        average_vertex_interaction_degree = vertex_interaction_degree_distribution.value_sequence.mean()
        beta_guess = 0.5 * average_vertex_interaction_degree * \
            (1. - gamma_guess) * (1. - gamma_prime_guess) / interaction_intensity

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.network_size = num_of_nodes
        self._parameters.interaction_intensity = interaction_intensity
        self._parameters.beta = beta_guess
        self._parameters.gamma = gamma_guess
        self._parameters.gamma_prime = gamma_prime_guess
        # pylint: enable=attribute-defined-outside-init

        print('\n'.join([
            '\nHypergraph model parameters after setting from data set:',
            f'network_size          = {self._parameters.network_size}',
            f'max_dim               = {self._parameters.max_dimension}',
            f'interaction_intensity = {self._parameters.interaction_intensity:4f}',
            f'beta                  = {self._parameters.beta:4f}',
            f'gamma                 = {self._parameters.gamma:4f}',
            f'gamma_prime           = {self._parameters.gamma_prime:4f}\n',
        ]))

    def generate_finite_network(self, seed: int) -> FiniteNetwork:
        """Build a network of the model."""
        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed}.')
        assert isinstance(self.parameters, HypergraphModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        self.parameters.torus_size_in_1_dimension = \
            self.parameters.network_size ** (1. / self.parameters.torus_dimension)

        if self.parameters.torus_dimension == 1:
            connections, interactions, nodes = generate_finite_network_cpp(self.parameters.to_numpy(), seed)
            node_birth_times, node_positions = nodes[:, 0], nodes[:, 1]
            interaction_birth_times, interaction_positions = interactions[:, 0], interactions[:, 1]
            interactions = pd.DataFrame(
                connections, columns=['interaction_id', 'vertex_id']
            ).groupby('interaction_id')['vertex_id'].apply(list).tolist()
        else:
            self._random_number_generator = np.random.default_rng(seed)
            node_birth_times, node_positions = self._create_vertex_points()
            interaction_birth_times, interaction_positions = self._create_interaction_points()
            interactions = self.generate_finite_network_interactions(
                node_birth_times, node_positions, interaction_birth_times, interaction_positions
            )

        vertex_ids = list(range(len(node_birth_times)))
        network = FiniteNetwork(self.parameters.max_dimension, vertex_ids, interactions)

        network.vertex_positions = self._create_vertex_positions_dict(node_positions, node_birth_times)
        network.interaction_positions = self._create_interaction_positions_dict(
            interaction_positions,
            interaction_birth_times,
        )

        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed} done.')
        return network

    def generate_infinite_network_set(self, num_of_networks: int, seed: int) -> InfiniteNetworkSet:

        if self.parameters.torus_dimension == 1:
            infinite_networks_cpp = generate_infinite_networks_cpp(
                self.parameters.to_numpy(),
                num_of_networks,
                seed
            )
            infinite_networks = []
            for (connections, interactions, vertices) in infinite_networks_cpp:
                interaction_birth_times, interaction_positions = interactions[:, 0], interactions[:, 1]
                vertex_birth_times, vertex_positions = vertices[:, 0], vertices[:, 1]

                vertex_ids = np.array(range(len(vertices)))[:, np.newaxis]
                interactions = pd.DataFrame(
                    connections, columns=['interaction_id', 'vertex_id']
                ).groupby('interaction_id')['vertex_id'].apply(list).tolist()
                network = InfiniteNetwork(self.parameters.max_dimension, vertex_ids, interactions)

                network.vertex_positions = self._create_vertex_positions_dict(
                    vertex_positions,
                    vertex_birth_times
                )
                network.interaction_positions = self._create_interaction_positions_dict(
                    interaction_positions,
                    interaction_birth_times,
                )

                infinite_networks.append(network)
        else:
            raise NotImplementedError

        return InfiniteNetworkSet(infinite_networks)

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
        for interaction_birth_time, interaction_position in tqdm(zip(
            interaction_birth_times.tolist(),
            interaction_positions.tolist(),
        ), delay=10, total=len(interaction_birth_times), desc='Generating interactions', leave=False):
            distances = self._distances(
                np.c_[np.full(len(node_birth_times), interaction_position), node_positions]
            )
            distances_d = distances**self.parameters.torus_dimension
            profile_function_argument = distances_d / self.parameters.beta \
                * node_birth_times**(self.parameters.gamma) \
                * interaction_birth_time**(self.parameters.gamma_prime)
            interaction_members = np.where(profile_function_argument < 1.)[0]
            interactions.append(interaction_members.tolist())

        return interactions

    def _create_vertex_points(self) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_]]:

        num_of_nodes = self._random_number_generator.poisson(lam=self.parameters.network_size)
        node_birth_times = self._random_number_generator.random(size=num_of_nodes)
        node_positions = self._random_number_generator.uniform(
            -self.parameters.torus_size_in_1_dimension / 2,
            +self.parameters.torus_size_in_1_dimension / 2,
            size=(num_of_nodes, self.parameters.torus_dimension)
        )
        return node_birth_times, node_positions

    def _create_interaction_points(self) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_]]:

        num_of_interactions = self._random_number_generator.poisson(
            lam=self.parameters.network_size * self.parameters.interaction_intensity)
        interaction_birth_times = self._random_number_generator.random(size=num_of_interactions)
        interaction_positions = self._random_number_generator.uniform(
            -self.parameters.torus_size_in_1_dimension / 2,
            +self.parameters.torus_size_in_1_dimension / 2,
            size=(num_of_interactions, self.parameters.torus_dimension)
        )
        return interaction_birth_times, interaction_positions

    def _create_vertex_positions_dict(
        self,
        node_positions: npt.NDArray[np.float64],
        node_birth_times: npt.NDArray[np.float64],
    ) -> dict[int, tuple[float, ...]]:
        node_ids = np.array(range(len(node_positions)))[:, np.newaxis]
        if self.parameters.torus_dimension == 1:
            vertex_positions_dict = {
                node_id: (node_position, node_birth_time)
                for node_id, node_position, node_birth_time
                in np.c_[node_ids, node_positions, node_birth_times]
            }
        elif self.parameters.torus_dimension == 2:
            vertex_positions_dict = {
                node_id: (node_position_0, node_position_1, node_birth_time)
                for node_id, node_position_0, node_position_1, node_birth_time
                in np.c_[node_ids, node_positions, node_birth_times]
            }
        return vertex_positions_dict

    def _create_interaction_positions_dict(
        self,
        interaction_positions: npt.NDArray[np.float64],
        interaction_birth_times: npt.NDArray[np.float64],
    ) -> dict[int, tuple[float, ...]]:
        interaction_ids = np.array(range(len(interaction_positions)))[:, np.newaxis]
        if self.parameters.torus_dimension == 1:
            interaction_positions_dict = {
                node_id: (node_position, node_birth_time)
                for node_id, node_position, node_birth_time
                in np.c_[interaction_ids, interaction_positions, interaction_birth_times]
            }
        elif self.parameters.torus_dimension == 2:
            interaction_positions_dict = {
                node_id: (node_position_0, node_position_1, node_birth_time)
                for node_id, node_position_0, node_position_1, node_birth_time
                in np.c_[interaction_ids, interaction_positions, interaction_birth_times]
            }
        return interaction_positions_dict

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
        self.parameters = HypergraphModel.Parameters(*parameters_tuple)

    def get_info_as_dict(self) -> dict[str, int | float]:
        """Return a dict representation based on the model properties."""
        return {
            'max_dimension': self.parameters.max_dimension,
            'network_size': self.parameters.network_size,
            'parameter_beta': self.parameters.beta,
            'parameter_gamma': self.parameters.gamma,
            'parameter_gamma_prime': self.parameters.gamma_prime,
            'torus_dimension': self.parameters.torus_dimension,
            'torus_size_in_1_dimension': self.parameters.torus_size_in_1_dimension,
        }

    @ property
    def parameters(self) -> HypergraphModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @ parameters.setter
    def parameters(self, value: HypergraphModel.Parameters) -> None:
        self._parameters = value
