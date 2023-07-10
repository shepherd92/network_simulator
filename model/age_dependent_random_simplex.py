#!/usr/bin/env python3
"""This module represents the age-dependent random simplex network model."""

from __future__ import annotations

from dataclasses import dataclass
from logging import info

import numpy as np
import numpy.typing as npt
from tqdm import tqdm

from cpp_modules.build.adrcm import (
    generate_finite_network_connections_default,
    generate_infinite_network_connections_default,
)
from data_set.data_set import DataSet
from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import create_fitting_parameters_power_law_data_set
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from model.model import Model
from network.finite_network import FiniteNetwork
from network.infinite_network import InfiniteNetwork, InfiniteNetworkSet
from network.property import BaseNetworkProperty


class AgeDependentRandomSimplexModel(Model):
    """Class representing an age-dependent random simplex network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an AgeDependentRandomSimplexModel."""

        torus_dimension: int = 0
        torus_size: float = 0.5  # the total size of the torus
        alpha: float = 0.5  # parameter of the profile function
        beta: float = 0.1  # edge density is beta / (1 - gamma)
        gamma: float = 0.5  # power law exponent is 1 + 1 / gamma

        def to_numpy(self) -> npt.NDArray[np.float_]:
            """Return the parameters as a numpy array."""
            return np.array([
                self.max_dimension,
                self.num_nodes,
                self.alpha,
                self.beta,
                self.gamma,
                self.torus_dimension,
                self.torus_size,
            ])

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = AgeDependentRandomSimplexModel.Parameters()

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
        fitting_parameters = create_fitting_parameters_power_law_data_set()
        approximation.fit(fitting_parameters)
        gamma_guess = 1. / (approximation.theoretical.parameters.exponent - 1.)
        beta_guess = (1. - gamma_guess) * average_degree

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.num_nodes = num_of_nodes
        # pylint: enable=attribute-defined-outside-init
        self._parameters.beta = beta_guess
        self._parameters.gamma = gamma_guess

        print('\n'.join([
            '\nADRCM model paramerers after setting from data set:',
            f'N       = {self._parameters.num_nodes}',
            f'max_dim = {self._parameters.max_dimension}',
            f'alpha   = {self._parameters.alpha:4f}',
            f'beta    = {self._parameters.beta:4f}',
            f'gamma   = {self._parameters.gamma:4f}\n',
        ]))

    def generate_finite_network(self, seed: int | None = None) -> FiniteNetwork:
        """Build a network of the model."""
        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed}.')
        assert isinstance(self.parameters, AgeDependentRandomSimplexModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        node_ids = np.array(range(self.parameters.num_nodes))[:, np.newaxis]
        random_number_generator = np.random.default_rng(seed)

        interarrival_times: npt.NDArray[np.float_] = \
            random_number_generator.exponential(size=self.parameters.num_nodes)
        time = interarrival_times.sum()
        birth_times = interarrival_times.cumsum() / time
        self.parameters.torus_size = time

        positions = random_number_generator.uniform(
            - time / 2, time / 2,
            size=(self.parameters.num_nodes, self.parameters.torus_dimension)
        )

        connections = self._generate_connections(birth_times, positions, seed)

        network = FiniteNetwork(self.parameters.max_dimension)
        network.add_simplices_batch(node_ids)
        network.add_simplices_batch(connections)
        network.expand()
        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed} done.')

        return network

    def generate_infinite_network_set(self, num_of_networks: int, seed: int) -> InfiniteNetworkSet:
        """Generate a list of infinite networks."""
        info(f'Generating infinite network set ({self.__class__.__name__}) with seed {seed}.')
        assert self.parameters.torus_dimension == 1, \
            f'Torus dimension must be 1, but it is {self.parameters.torus_dimension}'

        connections_set: list[npt.NDArray[np.float_]] = \
            generate_infinite_network_connections_default(
                self.parameters.to_numpy(),
                num_of_networks,
                seed
        )

        infinite_network_set = InfiniteNetworkSet([
            self.generate_infinite_network(connections)
            for connections in connections_set
        ])
        info(f'Generating infinite network set ({self.__class__.__name__}) with seed {seed} done.')
        return infinite_network_set

    def generate_infinite_network(self, connections: npt.NDArray[np.int_]) -> InfiniteNetwork:
        """Generate an "infinite" network, where the typical simplices are the ones that contain vertex 0."""
        # as each node that we consider is connected to O, the number of nodes is the maximal index of connections
        network = InfiniteNetwork(self.parameters.max_dimension)
        # network.digraph.add_nodes_from(node_ids)
        # network.digraph.add_edges_from(connections)
        # network.graph = network.digraph.to_undirected()
        network.add_simplex([0])  # ensure that the origin is added
        network.add_simplices_batch(connections)
        network.expand()

        return network

    def generate_infinite_network_bkp(self, seed: int | None = None) -> InfiniteNetwork:
        """Generate an "infinite" network, where the typical simplices are the ones that contain vertex 0."""
        assert self.parameters.torus_dimension == 1, \
            f'Torus dimension must be 1, but it is {self.parameters.torus_dimension}'
        random_number_generator = np.random.default_rng(seed)

        b = self.parameters.beta
        g = self.parameters.gamma

        u = random_number_generator.uniform(0., 1.)  # birth time of oldest node
        self.parameters.torus_size = b / u
        # Z = b / g * (u**(-g) - 1)  # expected number of points to be generated
        # N = random_number_generator.poisson(Z)  # number of points to be generted

        N = int(np.round(b * (1/u - 1)))  # total area of the rectangle
        birth_times = np.r_[np.array([u]), random_number_generator.uniform(u, 1., size=N)]
        positions = np.r_[np.array([0.]), random_number_generator.uniform(-b/u, +b/u, size=N)]

        # vertices closer to the origin than 1/2 * beta * u^(-gamma) * v^(gamma - 1) are connected to  the origin
        mask = np.abs(positions) < 0.5 * b * u**(-g) * birth_times**(g - 1)
        birth_times_connected_to_o = birth_times[mask]
        positions_connected_to_o = positions[mask]

        node_ids = np.array(range(len(birth_times_connected_to_o)))
        connections = self._generate_connections(birth_times_connected_to_o, positions_connected_to_o, seed)

        network = InfiniteNetwork(self.parameters.max_dimension)
        network.digraph.add_nodes_from(node_ids)
        network.digraph.add_edges_from(connections)
        network.graph = network.digraph.to_undirected()
        network.generate_clique_complex_from_graph()

        return network

    def set_model_parameters_from_tuple(self, parameters_tuple: tuple[int]) -> None:
        """Convert a tuple to ModelParamters. Used for model optimization."""
        self.parameters = AgeDependentRandomSimplexModel.Parameters(*parameters_tuple)

    def get_info_as_dict(self) -> dict[str, int | float]:
        """Return a dict representation based on the model properties."""
        return {
            'num_of_nodes': self.parameters.num_nodes,
            'max_dimension': self.parameters.max_dimension,
            'torus_dimension': self.parameters.torus_dimension,
            'torus_size': self.parameters.torus_size,
            'parameter_alpha': self.parameters.alpha,
            'parameter_beta': self.parameters.beta,
            'parameter_gamma': self.parameters.gamma,
        }

    def _generate_connections(
        self,
        birth_times: npt.NDArray[np.float_],
        positions: npt.NDArray[np.float_],
        seed: int | None
    ) -> npt.NDArray[np.float_]:
        is_default_profile_function = np.isclose(self.parameters.alpha, 0.5)
        is_default_connection_generation = is_default_profile_function and self.parameters.torus_dimension == 1
        if is_default_connection_generation:
            return self._generate_connections_cpp(birth_times, positions, seed)
        else:
            return self._generate_connections_python(birth_times, positions, seed)

    def _generate_connections_cpp(
        self,
        birth_times: npt.NDArray[np.float_],
        positions: npt.NDArray[np.float_],
        seed: int | None
    ) -> npt.NDArray[np.float_]:
        assert np.isclose(self.parameters.alpha, 0.5), \
            'In C++, only the default connection generation is implemented. ' + \
            f'Alpha parameter is not 0.5, but {self.parameters.alpha}.'
        assert self.parameters.torus_dimension == 1, \
            'In C++, only the default connection generation is implemented. ' + \
            f'Torus dimension is not 1, but {self.parameters.torus_dimension}.'

        connections = generate_finite_network_connections_default(
            birth_times,
            positions.flatten(),
            self.parameters.to_numpy()
        )
        return connections

    def _generate_connections_python(
        self,
        birth_times: npt.NDArray[np.float_],
        positions: npt.NDArray[np.float_],
        seed: int | None
    ) -> npt.NDArray[np.float_]:

        if self.parameters.num_nodes <= 1000:
            connections = self._generate_connections_vectorized(birth_times, positions, seed)
        else:
            connections = self._generate_connections_not_vectorized(birth_times, positions, seed)

        return connections

    def _generate_connections_vectorized(
        self,
        birth_times: npt.NDArray[np.float_],
        positions: npt.NDArray[np.float_],
        seed: int | None
    ) -> npt.NDArray[np.float_]:

        random_number_generator = np.random.default_rng(seed)

        distance_matrix = self._calc_torus_distance_matrix(positions)
        connection_probabilities = self._calc_connection_probabilities(birth_times, distance_matrix)
        random_numbers = random_number_generator.random(size=connection_probabilities.shape)
        adjacency_matrix = (random_numbers < connection_probabilities)
        connections = np.c_[np.where(adjacency_matrix)]

        return connections

    def _generate_connections_not_vectorized(
        self,
        birth_times: npt.NDArray[np.float_],
        positions: npt.NDArray[np.float_],
        seed
    ) -> npt.NDArray[np.float_]:
        """Calculate connections for large graphs (consumes lower memory)."""
        random_number_generator = np.random.default_rng(seed)
        node_ids = np.array(range(self.parameters.num_nodes))

        # condition on default vs non-default profile function for higher performance
        is_default_profile_function = np.isclose(self.parameters.alpha, 0.5)

        connection_list: list[np.ndarray] = []
        for source_node_id, source_birth_time, source_position in tqdm(zip(
            node_ids.tolist(), birth_times.tolist(), positions.tolist()
        ), desc='Generating connections:', delay=10, total=len(node_ids)):
            target_birth_times = birth_times[:source_node_id]
            birth_time_ratios = source_birth_time / target_birth_times

            target_positions = positions[:source_node_id]
            distances = self._distances(
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

    def _calc_torus_distance_matrix(self, positions: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:

        num_of_nodes = len(positions)

        # generate all possible pairs
        source_node_positions = np.repeat(positions, repeats=num_of_nodes, axis=0)
        target_node_positions = np.tile(positions, reps=(num_of_nodes, 1))
        all_position_pairs = np.c_[source_node_positions, target_node_positions]

        # calculate distances between all pairs
        distances = self._distances(all_position_pairs)
        distance_matrix = distances.reshape((num_of_nodes, num_of_nodes))

        return distance_matrix

    def _distances(self, position_pairs: npt.NDArray[np.float_]) -> float:

        # calculate the distances by dimension inside the boundaries without going around the edges
        distances_by_dimension_inside = np.abs(
            position_pairs[:, :self.parameters.torus_dimension] -
            position_pairs[:, self.parameters.torus_dimension:]
        )

        # if the distance in one dimension is greater than the half of the size of the torus,
        # then it is better to go around the edge, and the distance in that direction is
        # torus_size - distance_inside
        distances_by_dimension = np.where(
            distances_by_dimension_inside < 0.5 * self.parameters.torus_size,
            distances_by_dimension_inside,
            self.parameters.torus_size - distances_by_dimension_inside
        )

        # the final distances are the length of the difference vectors
        distances = np.linalg.norm(distances_by_dimension, axis=1)

        return distances

    def _calc_connection_probabilities(
        self,
        birth_times: npt.NDArray[np.float_],
        distance_matrix: npt.NDArray[np.float_]
    ) -> npt.NDArray[np.float_]:
        """Calculate the matrix of connection probabilities."""
        num_nodes = len(birth_times)
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

    def _profile_function(self, argument: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
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
