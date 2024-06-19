#!/usr/bin/env python3
"""Base class for network models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto, unique
from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm

from data_set.data_set import DataSet
from distribution.empirical_distribution import EmpiricalDistribution
from network.finite_network import FiniteNetwork
from network.infinite_network import InfiniteNetworkSet
from network.property import BaseNetworkProperty, DerivedNetworkProperty


class Model:
    """Base class representing a network model."""

    @unique
    class Type(Enum):
        """Specify the available model types."""

        ERDOS_RENYI = auto()
        WATTS_STROGATZ = auto()
        PREFERENTIAL_ATTACHMENT = auto()
        PRICE = auto()
        NETWORK_GEOMETRY_WITH_FLAVOR = auto()
        AGE_DEPENDENT_RANDOM_SIMPLEX = auto()
        AGE_DEPENDENT_RANDOM_HYPERGRAPH = auto()
        INVALID = auto()

    @unique
    class Mode(Enum):
        """Model testing mode."""
        FINITE: int = auto()
        INFINITE: int = auto()

    @dataclass
    class Parameters:
        """Base class representing the parameters."""

        max_dimension: int = 2
        network_size: float = 1000.

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        raise NotImplementedError

    def simulate(
        self,
        mode: Mode,
        scalar_property_params_to_calculate: list[DerivedNetworkProperty],
        num_of_simulations: int,
        initial_seed: int,
        num_of_infinite_networks_per_simulation: int = 0,
    ) -> dict[str, EmpiricalDistribution]:
        network_simulation_needed = any([
            not scalar_property_params.directly_calculated_from_model
            for scalar_property_params in scalar_property_params_to_calculate
        ])

        scalar_properties: dict[str, list[int | float]] = {
            property_type.name: []
            for property_type in scalar_property_params_to_calculate
        }

        for seed in tqdm(
            range(initial_seed, initial_seed + num_of_simulations),
            desc='Simulation'
        ):
            if network_simulation_needed:
                network = self.generate_finite_network(seed) \
                    if mode == Model.Mode.FINITE \
                    else self.generate_infinite_network_set(num_of_infinite_networks_per_simulation, seed)

            network_base_properties: dict[BaseNetworkProperty, Any] = {
                scalar_property_params.source_base_property: network.calc_base_property(scalar_property_params.source_base_property)
                for scalar_property_params in scalar_property_params_to_calculate
                if not scalar_property_params.directly_calculated_from_model
            }

            model_base_properties: dict[BaseNetworkProperty, Any] = {
                scalar_property_params.source_base_property: self.calc_base_property(
                    scalar_property_params.source_base_property, num_of_infinite_networks_per_simulation, seed)
                for scalar_property_params in scalar_property_params_to_calculate
                if scalar_property_params.directly_calculated_from_model
            }

            base_properties = {**network_base_properties, **model_base_properties}

            for scalar_property_params in scalar_property_params_to_calculate:
                scalar_properties[scalar_property_params.name].append(
                    scalar_property_params.calculator_default(
                        base_properties[scalar_property_params.source_base_property]
                    )
                )

        return {name: EmpiricalDistribution(values) for name, values in scalar_properties.items()}

    def generate_finite_network(self, seed: int) -> FiniteNetwork:
        """Build a network of the model."""
        raise NotImplementedError

    def generate_infinite_network_set(self, num_of_networks: int, seed: int) -> InfiniteNetworkSet:
        """Generate a set of "infinite" networks."""
        raise NotImplementedError

    def calc_base_property(
        self,
        property_type: BaseNetworkProperty,
        num_of_networks: int,
        seed: int
    ) -> Any:
        """Calculate the base property of the model."""
        raise NotImplementedError

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the model properties."""
        return {
            'network_size': self.parameters.network_size,
            'max_dimension': self.parameters.max_dimension,
        }

    def save_info(self, save_path: Path) -> None:
        """Save the main parameters to the given file as a pandas data frame."""
        information = self.info()
        data_frame = pd.DataFrame(information, index=[0])
        data_frame.to_csv(save_path, index=False)

    @property
    def parameters(self) -> Model.Parameters:
        """Get the simple graph associated to the data set."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: Model.Parameters) -> None:
        self._parameters = value

    @parameters.deleter
    def parameters(self) -> None:
        del self._parameters

    @property
    def parameters_type(self) -> type:
        """Return the model parameters type."""
        return self.Parameters

    def __str__(self) -> str:
        """Return a string representation based on the network properties."""
        return '\n'.join([
            f'{key}: {item}'
            for key, item in self.info().items()
        ])
