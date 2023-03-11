#!/usr/bin/env python3
"""This module represents the age-dependent random simplex network model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from data_set.data_set import DataSet
from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from model.model import Model
from network.network import Network
from network.property import BaseNetworkProperty


class AgeDependentRandomSimplexModel(Model):
    """Class representing an age-dependent random simplex network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an AgeDependentRandomSimplexModel."""

        torus_dimension: int = 0
        alpha: float = 0.5  # parameter of the profile function
        beta: float = 0.1  # edge density is beta / (1 - gamma)
        gamma: float = 0.5  # power law exponent is 1 + 1 / gamma

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters: AgeDependentRandomSimplexModel.Parameters = AgeDependentRandomSimplexModel.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)
        average_degree: float = data_set.calc_base_property(BaseNetworkProperty.Type.AVERAGE_DEGREE)
        degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION
        )
        approximation = DistributionApproximation(
            degree_distribution,
            TheoreticalDistribution.Type.POWER_LAW
        )
        approximation.fit()
        gamma_guess = 1. / (approximation.theoretical.parameters.exponent - 1.)
        beta_guess = (1. - gamma_guess) * average_degree

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.num_nodes = num_of_nodes
        # pylint: enable=attribute-defined-outside-init
        self._parameters.beta = beta_guess
        self._parameters.gamma = gamma_guess

    def generate_network(self, seed: int | None = None) -> Network:
        """Build a network of the model."""
        assert isinstance(self.parameters, AgeDependentRandomSimplexModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        node_ids = np.array(range(self.parameters.num_nodes))
        random_number_generator = np.random.default_rng(seed)

        interarrival_times: np.ndarray = random_number_generator.exponential(size=self.parameters.num_nodes)
        birth_times = interarrival_times.cumsum()
        positions = random_number_generator.random(
            size=(self.parameters.num_nodes, self.parameters.torus_dimension)
        )

        if self.parameters.num_nodes <= 1000:
            connections = self._generate_connections_vectorized(birth_times, positions, seed)
        else:
            connections = self._generate_connections_not_vectorized(birth_times, positions, seed)

        network = Network(self.parameters.max_dimension)
        network.digraph.add_nodes_from(node_ids)
        network.digraph.add_edges_from(connections)
        network.graph = network.digraph.to_undirected()
        network.generate_clique_complex_from_graph()

        return network

    def set_model_parameters_from_tuple(self, parameters_tuple: tuple[int]) -> None:
        """Convert a tuple to ModelParamters. Used for model optimization."""
        self.parameters = AgeDependentRandomSimplexModel.Parameters(*parameters_tuple)

    def _generate_connections_vectorized(
        self,
        birth_times: np.ndarray,
        positions: np.ndarray,
        seed: int | None
    ) -> np.ndarray:

        random_number_generator = np.random.default_rng(seed)

        distance_matrix = AgeDependentRandomSimplexModel._calc_torus_distance_matrix(positions)
        connection_probabilities = self._calc_connection_probabilities(birth_times, distance_matrix)
        random_numbers = random_number_generator.random(size=connection_probabilities.shape)
        adjacency_matrix = (random_numbers < connection_probabilities)
        connections = np.c_[np.where(adjacency_matrix)]

        return connections

    def _generate_connections_not_vectorized(
        self,
        birth_times: np.ndarray,
        positions: np.ndarray,
        seed
    ) -> np.ndarray:
        """Calculate connections for large graphs (consumes lower memory)."""
        random_number_generator = np.random.default_rng(seed)
        node_ids = np.array(range(self.parameters.num_nodes))

        # condition on default vs non-default profile function for higher performance
        is_default_profile_function = np.isclose(self.parameters.alpha, 0.5)

        connection_list: list[np.ndarray] = []
        for source_node_id, source_birth_time, source_position in zip(
            node_ids.tolist(), birth_times.tolist(), positions.tolist()
        ):
            target_birth_times = birth_times[:source_node_id]
            birth_time_ratios = source_birth_time / target_birth_times

            target_positions = positions[:source_node_id]
            distances = AgeDependentRandomSimplexModel._distances(
                np.c_[np.full(source_node_id, source_position), target_positions]
            )
            distances_d_x_birth_times = distances**self.parameters.torus_dimension * source_birth_time
            profile_function_argument = distances_d_x_birth_times / \
                (self.parameters.beta * birth_time_ratios**self.parameters.gamma)

            if is_default_profile_function:
                # use the default profile function, deterministic
                adjacency_matrix_row = profile_function_argument <= 0.5
            else:
                # use the general profile function
                connection_probabilities = self._profile_function(profile_function_argument)
                adjacency_matrix_row = (
                    random_number_generator.random(size=connection_probabilities.shape) <
                    connection_probabilities
                )
            target_node_ids: npt.NDArray[np.float_] = np.c_[np.where(adjacency_matrix_row)]
            connections_from_source_node = np.c_[
                np.full((target_node_ids.shape), source_node_id),
                target_node_ids
            ]
            connection_list.append(connections_from_source_node)

        connections = np.concatenate(connection_list)

        return connections

    def _generate_typical_in_degree(self) -> int:
        beta = self.parameters.beta
        gamma = self.parameters.gamma
        origin_birth_time = np.random.uniform()
        network_size = 1. / origin_birth_time
        birth_times = np.random.uniform(0., 1., int(network_size))
        birth_times_of_younger_nodes = birth_times[birth_times > origin_birth_time]
        num_of_younger_nodes = len(birth_times_of_younger_nodes)
        positions = np.random.uniform(-network_size / 2., +network_size / 2., num_of_younger_nodes)

        profile_function_arguments = birth_times_of_younger_nodes * np.abs(positions) / \
            (beta * (birth_times_of_younger_nodes / origin_birth_time)**gamma)
        in_connections = self._profile_function(profile_function_arguments)
        in_degree = in_connections.sum()

        return in_degree

    @staticmethod
    def _calc_torus_distance_matrix(positions: np.ndarray) -> np.ndarray:

        num_of_nodes = len(positions)

        # generate all possible pairs
        source_node_positions = np.repeat(positions, repeats=num_of_nodes, axis=0)
        target_node_positions = np.tile(positions, reps=(num_of_nodes, 1))
        all_position_pairs = np.c_[source_node_positions, target_node_positions]

        # calculate distances between all pairs
        distances = AgeDependentRandomSimplexModel._distances(all_position_pairs)
        distance_matrix = distances.reshape((num_of_nodes, num_of_nodes))

        return distance_matrix

    @staticmethod
    def _distances(position_pairs: np.ndarray) -> float:

        dimensions = position_pairs.shape[1] // 2

        # calculate the distances by dimension inside the boundaries without going around the edges
        distances_by_dimension_inside = np.abs(position_pairs[:, :dimensions] - position_pairs[:, dimensions:])

        # if the distance in one dimension is greater than the half of the size of the torus,
        # then it is better to go around the edge, and the distance in that direction is
        # 1 - distance_inside
        distances_by_dimension = np.where(
            distances_by_dimension_inside < 0.5,
            distances_by_dimension_inside,
            1 - distances_by_dimension_inside
        )

        # the final distances are the length of the difference vectors
        distances = np.linalg.norm(distances_by_dimension, axis=1)

        return distances

    def _calc_connection_probabilities(
        self,
        birth_times: np.ndarray,
        distance_matrix: np.ndarray
    ) -> np.ndarray:
        """Calculate the matrix of connection probabilities."""
        num_nodes = self.parameters.num_nodes
        source_node_birth_times = np.repeat(birth_times, repeats=num_nodes, axis=0)
        target_node_birth_times = np.tile(birth_times, reps=(num_nodes))
        all_birth_time_pairs = np.c_[source_node_birth_times, target_node_birth_times]
        birth_time_ratios = all_birth_time_pairs[:, 0] / all_birth_time_pairs[:, 1]
        birth_time_ratio_matrix = birth_time_ratios.reshape((num_nodes, num_nodes))

        distances_d_x_birth_times = distance_matrix**self.parameters.torus_dimension * birth_times[:, np.newaxis]
        profile_function_argument = distances_d_x_birth_times / \
            (self.parameters.beta * birth_time_ratio_matrix**self.parameters.gamma)

        connection_probabilities = self._profile_function(profile_function_argument)

        # source node birth time must be greater than target node birth time:
        # source nodes are born after target nodes, these entries are in the lower triangle matrix
        # the only valid probabilities are thus in the lower triangle
        connection_probabilities *= np.tri(*connection_probabilities.shape, k=-1)

        return connection_probabilities

    def _profile_function(self, argument: np.ndarray) -> np.ndarray:
        result = np.zeros_like(argument)
        result[argument <= self.parameters.alpha] = 1. / (2. * self.parameters.alpha)
        return result

    @property
    def parameters(self) -> AgeDependentRandomSimplexModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: AgeDependentRandomSimplexModel.Parameters) -> None:
        self._parameters = value
