#!/usr/bin/env python3
"""This module represents the age-dependent random hypergraph network model."""

from __future__ import annotations

from dataclasses import dataclass
from logging import info

import numpy as np
import numpy.typing as npt
from scipy import integrate
from scipy.optimize import fsolve

# pylint: disable-next=no-name-in-module
from cpp_plugin.build.release.cpp_plugin import (  # type: ignore
    FiniteHypergraphModel,
    InfiniteHypergraphModel,
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

        def to_numpy(self) -> npt.NDArray[np.float_]:
            """Return the parameters as a numpy array."""
            return np.array([
                self.max_dimension,
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
        num_of_vertices: int = data_set.calc_base_property(BaseNetworkProperty.num_of_vertices)
        num_of_interactions: int = data_set.calc_base_property(BaseNetworkProperty.num_of_interactions)

        vertex_interaction_degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.vertex_interaction_degree_distribution
        )
        vertex_interaction_exponent = guess_power_law_exponent(vertex_interaction_degree_distribution)
        average_vertex_interaction_degree = vertex_interaction_degree_distribution.value_sequence.mean()
        gamma = 1. / (vertex_interaction_exponent - 1.)
        assert gamma > 0. and gamma < 1., f'gamma_guess = {gamma}'

        document_degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.interaction_vertex_degree_distribution
        )
        interaction_vertex_exponent = guess_power_law_exponent(document_degree_distribution)
        gamma_prime = 1. / (interaction_vertex_exponent - 1.)
        assert gamma_prime > 0. and gamma_prime < 1., f'gamma_prime_guess = {gamma_prime}'

        def system_of_equations(x: npt.NDArray[np.float_]) -> float:
            network_size = x[0]
            interaction_intensity = x[1]

            beta = 0.5 * average_vertex_interaction_degree * (1. - gamma) * (1. - gamma_prime) / interaction_intensity

            def integrand(y: float, lambda_prime: float, g: float, g_prime: float) -> float:
                return np.exp(-2. * beta * lambda_prime / (1. - g_prime) * y**(-g))
            vertex_correction_integral = integrate.quad(
                integrand, 0., 1., args=(interaction_intensity, gamma, gamma_prime))[0]
            interaction_correction_integral = integrate.quad(
                integrand, 0., 1., args=(network_size, gamma_prime, gamma))[0]

            num_vertices_error = num_of_vertices - network_size * (1. - vertex_correction_integral)
            num_interactions_error = \
                num_of_interactions - interaction_intensity * (1. - interaction_correction_integral)

            return [num_vertices_error, num_interactions_error]

        network_size, interaction_intensity = fsolve(
            system_of_equations, (num_of_vertices, num_of_interactions),
        )
        beta_guess = 0.5 * average_vertex_interaction_degree * (1. - gamma) * (1. - gamma_prime) / interaction_intensity

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.network_size = network_size
        self._parameters.interaction_intensity = interaction_intensity
        self._parameters.beta = beta_guess
        self._parameters.gamma = gamma
        self._parameters.gamma_prime = gamma_prime
        # pylint: enable=attribute-defined-outside-init

        print(self)

    def generate_finite_network(self, seed: int) -> FiniteNetwork:
        """Build a network of the model."""
        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed}.')

        cpp_model = FiniteHypergraphModel(self._parameters.to_numpy(), seed)
        cpp_network, vertex_positions, interaction_positions = cpp_model.generate_network()
        network = FiniteNetwork(cpp_network)
        network.vertex_positions = self._create_point_positions_dict(vertex_positions)
        network.interaction_positions = self._create_point_positions_dict(interaction_positions)

        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed} done.')
        return network

    def generate_infinite_network_set(self, num_of_networks: int, seed: int) -> InfiniteNetworkSet:
        """Build infinite networks of the model."""
        info(f'Generating infinite network set ({self.__class__.__name__}) with seed {seed}.')

        cpp_model = InfiniteHypergraphModel(self._parameters.to_numpy(), seed)
        cpp_networks = cpp_model.generate_networks(num_of_networks)

        infinite_networks: list[InfiniteNetwork] = []
        for cpp_network, vertex_positions, interaction_positions in cpp_networks:
            infinite_network = InfiniteNetwork(cpp_network)
            infinite_network.vertex_positions = self._create_point_positions_dict(vertex_positions)
            infinite_network.interaction_positions = self._create_point_positions_dict(interaction_positions)
            infinite_networks.append(InfiniteNetwork(cpp_network))
        return InfiniteNetworkSet(infinite_networks)

    def _create_point_positions_dict(self, vertex_positions: list[tuple[float, float]]) -> dict[int, tuple[float, ...]]:
        vertex_ids = list(range(len(vertex_positions)))
        vertex_positions_dict = {
            id: (position, mark)
            for id, (mark, position)
            in zip(vertex_ids, vertex_positions)
        }
        return vertex_positions_dict

    def set_model_parameters_from_tuple(self, parameters_tuple: tuple[int]) -> None:
        """Convert a tuple to ModelParamters. Used for model optimization."""
        self.parameters = HypergraphModel.Parameters(*parameters_tuple)

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the model properties."""
        return {
            'max_dimension': self.parameters.max_dimension,
            'network_size': self.parameters.network_size,
            'interaction_intensity': self.parameters.interaction_intensity,
            'parameter_beta': self.parameters.beta,
            'parameter_gamma': self.parameters.gamma,
            'parameter_gamma_prime': self.parameters.gamma_prime,
        }

    @ property
    def parameters(self) -> HypergraphModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @ parameters.setter
    def parameters(self, value: HypergraphModel.Parameters) -> None:
        self._parameters = value
