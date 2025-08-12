#!/usr/bin/env python3
"""This module represents the age-dependent random simplex network model."""

from __future__ import annotations

from dataclasses import dataclass
from logging import info

import numpy as np
import numpy.typing as npt

# pylint: disable-next=no-name-in-module
from cpp_plugin.build.release.cpp_plugin import (
    FiniteAdrcmModel,
    InfiniteAdrcmModel,
)
from data_set.data_set import DataSet
from distribution.approximation import guess_power_law_exponent
from distribution.empirical_distribution import EmpiricalDistribution
from model.model import Model
from network.finite_clique_complex import FiniteCliqueComplex
from network.infinite_clique_complex import InfiniteCliqueComplex, InfiniteCliqueComplexSet
from network.property import BaseNetworkProperty


class AdrcmModel(Model):
    """Class representing an age-dependent random simplex network model."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an AdrcmModel."""

        alpha: float = 0.5  # parameter of the profile function
        beta: float = 0.1  # edge density is beta / (1 - gamma)
        gamma: float = 0.5  # power law exponent is 1 + 1 / gamma

        def to_numpy(self) -> npt.NDArray[np.float64]:
            """Return the parameters as a numpy array."""
            return np.array([
                self.max_dimension,
                self.network_size,
                self.alpha,
                self.beta,
                self.gamma,
            ])

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        super().__init__()
        self._parameters = AdrcmModel.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.num_of_vertices)
        average_degree: float = data_set.calc_base_property(BaseNetworkProperty.mean_degree)
        degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.vertex_edge_degree_distribution
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

        print(self)

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the model properties."""
        return {
            'max_dimension': self.parameters.max_dimension,
            'network_size': self.parameters.network_size,
            'parameter_alpha': self.parameters.alpha,
            'parameter_beta': self.parameters.beta,
            'parameter_gamma': self.parameters.gamma,
        }

    def generate_finite_network(self, seed: int) -> FiniteCliqueComplex:
        """Build a network of the model."""
        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed}.')
        assert isinstance(self.parameters, AdrcmModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'

        cpp_model = FiniteAdrcmModel(self._parameters.to_numpy(), seed)
        cpp_network, vertex_positions = cpp_model.generate_network()
        network = FiniteCliqueComplex(cpp_network)
        network.vertex_positions = self._create_point_positions_dict(vertex_positions)

        info(f'Generating finite network ({self.__class__.__name__}) with seed {seed} done.')
        return network

    def generate_infinite_network_set(self, num_of_networks: int, seed: int) -> InfiniteCliqueComplexSet:
        """Generate a list of infinite networks."""
        info(f'Generating infinite network set ({self.__class__.__name__}) with seed {seed}.')

        cpp_model = InfiniteAdrcmModel(self._parameters.to_numpy(), seed)
        networks: list[InfiniteCliqueComplex] = [
            InfiniteCliqueComplex(cpp_network)
            for cpp_network in cpp_model.generate_networks(num_of_networks)
        ]
        infinite_network_set = InfiniteCliqueComplexSet(networks)

        info(f'Generating infinite network set ({self.__class__.__name__}) with seed {seed} done.')
        return infinite_network_set

    def set_model_parameters_from_tuple(self, parameters_tuple: tuple[int]) -> None:
        """Convert a tuple to Model.Parameters. Used for model optimization."""
        self.parameters = AdrcmModel.Parameters(*parameters_tuple)

    def _create_point_positions_dict(self, vertex_positions: list[tuple[float, float]]) -> dict[int, tuple[float, ...]]:
        vertex_ids = list(range(len(vertex_positions)))
        vertex_positions_dict = {
            id: (position, mark)
            for id, (mark, position)
            in zip(vertex_ids, vertex_positions)
        }
        return vertex_positions_dict

    @property
    def parameters(self) -> AdrcmModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: AdrcmModel.Parameters) -> None:
        self._parameters = value
