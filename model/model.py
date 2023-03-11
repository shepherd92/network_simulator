#!/usr/bin/env python3
"""Base class for network models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from multiprocessing import Pool
from typing import Any

from tqdm import tqdm

from data_set.data_set import DataSet
from distribution.empirical_distribution import EmpiricalDistribution
from network.network import Network
from network.property import BaseNetworkProperty, DerivedNetworkProperty
import tools.istarmap  # pylint: disable=unused-import # noqa: F401


class Model:
    """Base class representing a network model."""

    class Type(Enum):
        """Specify the available model types."""

        ERDOS_RENYI = auto()
        WATTS_STROGATZ = auto()
        PREFERENTIAL_ATTACHMENT = auto()
        PRICE = auto()
        AGE_DEPENDENT_RANDOM_CONNECTION = auto()
        NETWORK_GEOMETRY_WITH_FLAVOR = auto()
        AGE_DEPENDENT_RANDOM_SIMPLEX = auto()
        INVALID = auto()

    @dataclass
    class Parameters:
        """Base class representing the parameters."""

        max_dimension: int = 0
        num_nodes: int = 0

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = Model.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        raise NotImplementedError

    def simulate(
        self,
        scalar_property_params_to_calculate: list[DerivedNetworkProperty],
        num_of_simulations: int,
        num_of_processes: int
    ) -> list[EmpiricalDistribution]:
        """Calculate the summaries for the given network model."""
        base_network_properties = [
            scalar_property_params.source_base_property
            for scalar_property_params in scalar_property_params_to_calculate
        ]

        all_networks_base_network_properties = self._simulate_base_properties(
            base_network_properties,
            num_of_simulations,
            num_of_processes
        )

        scalar_property_distributions = self._extract_derived_properties(
            scalar_property_params_to_calculate,
            all_networks_base_network_properties
        )
        return scalar_property_distributions

    def _simulate_base_properties(
        self,
        base_network_properties: list[BaseNetworkProperty],
        num_of_simulations: int,
        num_of_processes: int
    ) -> list[list[Any]]:
        if num_of_processes > 1:
            base_network_property_values: list[list[Any]] = self._simulate_multiple_processes(
                base_network_properties,
                num_of_simulations,
                num_of_processes
            )
        else:
            base_network_property_values = self._simulate_single_process(
                base_network_properties,
                num_of_simulations
            )

        return base_network_property_values

    def _extract_derived_properties(
        self,
        scalar_property_params_to_calculate: list[DerivedNetworkProperty],
        base_network_property_values: list[list[Any]]
    ) -> list[EmpiricalDistribution]:
        scalar_property_distributions: list[EmpiricalDistribution] = []
        for index, scalar_property_params in enumerate(scalar_property_params_to_calculate):
            scalar_property_values = [
                scalar_property_params.calculator(this_network_base_properties[index])
                for this_network_base_properties in base_network_property_values
            ]
            empirical_distribution = EmpiricalDistribution(scalar_property_values)
            scalar_property_distributions.append(empirical_distribution)

        return scalar_property_distributions

    def generate_properties(self, base_properties: list[BaseNetworkProperty], seed: int) -> list[Any]:
        """Build a single network of the model and return its summary."""
        network: Network | None = None
        property_values: list[Any] = []
        for property_ in base_properties:
            if property_.calculation_method == BaseNetworkProperty.CalculationMethod.NETWORK:
                network = self.generate_network(seed) if network is None else network
                property_value = network.calc_base_property(property_.property_type)
            elif property_.calculation_method == BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT:
                property_value = self._calc_typical_property_distribution(property_.property_type)
            else:
                raise NotImplementedError(f'Unknown calculation method {property_.calculation_method}')
            property_values.append(property_value)

        return property_values

    def generate_network(self, seed: int | None = None) -> Network:
        """Build a network of the model."""
        raise NotImplementedError

    def _calc_typical_property_distribution(self, property_type: BaseNetworkProperty.Type) -> EmpiricalDistribution:
        """Generate typical properties of the given type."""
        num_of_objects = self.parameters.num_nodes
        generated_values: list[float | int] = [
            self._generate_typical_property(property_type)
            for _ in range(num_of_objects)
        ]

        distribution = EmpiricalDistribution(generated_values)
        return distribution

    def _simulate_multiple_processes(
        self,
        base_network_properties: list[BaseNetworkProperty],
        num_of_simulations: int,
        num_of_processes: int
    ) -> list[list[Any]]:
        """Simulate networks using multiple processes calculating a list of base properties."""
        with Pool(num_of_processes) as pool:
            # pylint: disable-next=no-member
            all_networks_base_network_properties: list[list[Any]] = tqdm(pool.istarmap(  # type: ignore
                self.generate_properties,
                list(zip([base_network_properties] * num_of_simulations, range(num_of_simulations)))
            ), total=num_of_simulations)
        return all_networks_base_network_properties

    def _simulate_single_process(
        self,
        base_network_properties: list[BaseNetworkProperty],
        num_of_simulations: int
    ) -> list[list[Any]]:
        all_networks_base_network_properties = [
            self.generate_properties(base_network_properties, seed)
            for seed in tqdm(range(num_of_simulations), total=num_of_simulations)
        ]
        return all_networks_base_network_properties

    def _generate_typical_property(self, property_type: BaseNetworkProperty.Type) -> float | int:
        if property_type == BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION:
            return self._generate_typical_in_degree()
        raise NotImplementedError(f'Unknown typical property {property_type}.')

    def _generate_typical_in_degree(self) -> int:
        raise NotImplementedError

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
