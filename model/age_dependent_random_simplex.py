#!/usr/bin/env python3
"""This module represents the age-dependent random simplex network model."""

from __future__ import annotations

from dataclasses import dataclass
from logging import info

import numpy as np
import numpy.typing as npt

from cpp_modules.build.adrcm import (
    generate_finite_network_connections,
    generate_infinite_network_connections,
)
from data_set.data_set import DataSet
from distribution.approximation import guess_power_law_exponent
from distribution.empirical_distribution import EmpiricalDistribution
from model.model import Model
from network.finite_network import FiniteNetwork
from network.infinite_network import InfiniteNetwork, InfiniteNetworkSet
from network.property import BaseNetworkProperty


class AgeDependentRandomSimplexModel(Model):
    """Class representing an age-dependent random simplex network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an AgeDependentRandomSimplexModel."""

        num_of_nodes: int = 0  # actual number of nodes
        alpha: float = 0.5  # parameter of the profile function
        beta: float = 0.1  # edge density is beta / (1 - gamma)
        gamma: float = 0.5  # power law exponent is 1 + 1 / gamma
        torus_dimension: int = 0
        torus_size_in_1_dimension: float = 0.5  # the total size of the torus

        def to_numpy(self) -> npt.NDArray[np.float_]:
            """Return the parameters as a numpy array."""
            return np.array([
                self.max_dimension,
                self.network_size,
                self.num_of_nodes,
                self.alpha,
                self.beta,
                self.gamma,
                self.torus_dimension,
                self.torus_size_in_1_dimension,
            ])

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        super().__init__()
        self._parameters = AgeDependentRandomSimplexModel.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)
        average_degree: float = data_set.calc_base_property(BaseNetworkProperty.Type.AVERAGE_DEGREE)
        degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION
        )
        exponent_guess = guess_power_law_exponent(degree_distribution)
        gamma_guess = 1. / (exponent_guess - 1.)
        beta_guess = (1. - gamma_guess) * average_degree

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.network_size = num_of_nodes
        # pylint: enable=attribute-defined-outside-init
        self._parameters.beta = beta_guess
        self._parameters.gamma = gamma_guess

        print('\n'.join([
            '\nADRCM model paramerers after setting from data set:',
            f'size    = {self._parameters.network_size}',
            f'max_dim = {self._parameters.max_dimension}',
            f'alpha   = {self._parameters.alpha:4f}',
            f'beta    = {self._parameters.beta:4f}',
            f'gamma   = {self._parameters.gamma:4f}\n',
        ]))

    def generate_finite_network(self, seed: int) -> FiniteNetwork:
        """Build a network of the model."""
        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed}.')
        assert isinstance(self.parameters, AgeDependentRandomSimplexModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        random_number_generator = np.random.default_rng(seed)
        self.parameters.num_of_nodes = random_number_generator.poisson(lam=self.parameters.network_size)
        self.parameters.torus_size_in_1_dimension = \
            self.parameters.network_size ** (1. / self.parameters.torus_dimension)

        node_ids = np.array(range(self.parameters.num_of_nodes))[:, np.newaxis]
        connections = generate_finite_network_connections(
            self.parameters.to_numpy(),
            seed,
        )

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
            generate_infinite_network_connections(
                self.parameters.to_numpy(),
                num_of_networks,
                seed,
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
        network.add_simplices([[0]])  # ensure that the origin is added
        network.add_simplices_batch(connections)
        network.expand()

        return network

    def set_model_parameters_from_tuple(self, parameters_tuple: tuple[int]) -> None:
        """Convert a tuple to ModelParamters. Used for model optimization."""
        self.parameters = AgeDependentRandomSimplexModel.Parameters(*parameters_tuple)

    def get_info_as_dict(self) -> dict[str, int | float]:
        """Return a dict representation based on the model properties."""
        return {
            'max_dimension': self.parameters.max_dimension,
            'network_size': self.parameters.network_size,
            'num_of_nodes': self.parameters.num_of_nodes,
            'parameter_alpha': self.parameters.alpha,
            'parameter_beta': self.parameters.beta,
            'parameter_gamma': self.parameters.gamma,
            'torus_dimension': self.parameters.torus_dimension,
            'torus_size_in_1_dimension': self.parameters.torus_size_in_1_dimension,
        }

    @property
    def parameters(self) -> AgeDependentRandomSimplexModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: AgeDependentRandomSimplexModel.Parameters) -> None:
        self._parameters = value
