#!/usr/bin/env python3
"""Base class for network models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from logging import info
from multiprocessing import Pool, Value
from tqdm import tqdm
from typing import Any

from data_set.data_set import DataSet
from distribution.empirical_distribution import EmpiricalDistribution
from network.finite_network import FiniteNetwork
from network.infinite_network import InfiniteNetwork
from network.property import BaseNetworkProperty, DerivedNetworkProperty
import tools.istarmap  # pylint: disable=unused-import # noqa: F401


NUM_OF_SIMULATIONS = 0
SIMULATIONS_DONE = Value('i', 0)
PROGRESS_BAR = None


class Model:
    """Base class representing a network model."""

    class Type(Enum):
        """Specify the available model types."""

        ERDOS_RENYI = auto()
        WATTS_STROGATZ = auto()
        PREFERENTIAL_ATTACHMENT = auto()
        PRICE = auto()
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
        num_of_infinite_networks: int,
        num_of_processes: int,
    ) -> list[EmpiricalDistribution]:
        """Calculate the summaries for the given network model."""
        base_network_properties = [
            scalar_property_params.source_base_property
            for scalar_property_params in scalar_property_params_to_calculate
        ]

        all_networks_base_network_properties = self._simulate_base_properties(
            base_network_properties,
            num_of_simulations,
            num_of_infinite_networks,
            num_of_processes,
        )

        scalar_property_distributions = self._extract_derived_properties(
            scalar_property_params_to_calculate,
            all_networks_base_network_properties
        )
        return scalar_property_distributions

    def generate_finite_network(self, seed: int) -> FiniteNetwork:
        """Build a network of the model."""
        raise NotImplementedError

    def generate_infinite_network_set(self, num_of_networks: int, seed: int) -> list[InfiniteNetwork]:
        """Generate a set of "infinite" networks."""
        raise NotImplementedError

    def generate_infinite_network(self, seed: int | None) -> InfiniteNetwork:
        """Generate an "infinite" network, where the typical simplices are the ones that contain vertex 0."""
        raise NotImplementedError

    def _simulate_base_properties(
        self,
        base_network_properties: list[BaseNetworkProperty],
        num_of_simulations: int,
        num_of_infinite_networks: int,
        num_of_processes: int
    ) -> list[list[Any]]:

        global NUM_OF_SIMULATIONS
        NUM_OF_SIMULATIONS = num_of_simulations

        global PROGRESS_BAR
        PROGRESS_BAR = tqdm(total=NUM_OF_SIMULATIONS, desc='Simulation')

        if num_of_processes > 1:
            base_network_property_values: list[list[Any]] = self._simulate_multiple_processes(
                base_network_properties,
                num_of_simulations,
                num_of_infinite_networks,
                num_of_processes,
            )
        else:
            base_network_property_values = self._simulate_single_process(
                base_network_properties,
                num_of_simulations,
                num_of_infinite_networks,
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

    def _simulate_multiple_processes(
        self,
        base_network_properties: list[BaseNetworkProperty],
        num_of_simulations: int,
        num_of_infinite_networks: int,
        num_of_processes: int,
    ) -> list[list[Any]]:
        """Simulate networks using multiple processes calculating a list of base properties."""
        info(f'Starting a pool of {num_of_processes} processes.')

        with Pool(num_of_processes) as pool:
            args = list(zip(
                [base_network_properties] * num_of_simulations,
                [num_of_infinite_networks] * num_of_simulations,
                range(num_of_simulations),
            ))
            # pylint: disable-next=no-member
            # all_networks_base_network_properties: list[list[Any]] = [
            #     results for results in
            #     tqdm(pool.istarmap(self._generate_properties, args), desc='Simulation:', total=num_of_simulations)
            # ]

            all_networks_base_network_properties: list[list[Any]] = pool.starmap(  # type: ignore
                self._generate_properties,
                args
            )
        info(f'Multiprocessing with {num_of_processes} processes finished.')
        return all_networks_base_network_properties

    def _simulate_single_process(
        self,
        base_network_properties: list[BaseNetworkProperty],
        num_of_simulations: int,
        num_of_infinite_networks: int,
    ) -> list[list[Any]]:
        all_networks_base_network_properties = [
            self._generate_properties(base_network_properties, num_of_infinite_networks, seed)
            for seed in range(num_of_simulations)
        ]
        return all_networks_base_network_properties

    def _generate_properties(
        self,
        base_properties: list[BaseNetworkProperty],
        num_of_infinite_networks: int,
        seed: int,
    ) -> list[Any]:
        """Build a single network of the model and return its summary."""
        finite_network: FiniteNetwork | None = None
        infinite_networks: list[InfiniteNetwork] | None = None
        property_values: list[Any] = []
        for property_ in base_properties:
            if property_.calculation_method == BaseNetworkProperty.CalculationMethod.NETWORK:
                finite_network = self.generate_finite_network(seed) if finite_network is None else finite_network
                property_value = finite_network.calc_base_property(property_.property_type)
            elif property_.calculation_method == BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT:
                infinite_networks = \
                    self.generate_infinite_network_set(num_of_infinite_networks, seed) \
                    if infinite_networks is None else infinite_networks
                property_value = self._calc_typical_property_distribution(infinite_networks, property_.property_type)
            else:
                raise NotImplementedError(f'Unknown calculation method {property_.calculation_method}')
            property_values.append(property_value)

        global SIMULATIONS_DONE
        global PROGRESS_BAR
        with SIMULATIONS_DONE.get_lock():
            SIMULATIONS_DONE.value += 1
            PROGRESS_BAR.n = SIMULATIONS_DONE.value
            PROGRESS_BAR.refresh()

        return property_values

    def _calc_typical_property_distribution(
        self,
        infinite_networks: list[InfiniteNetwork],
        property_type: BaseNetworkProperty.Type,
    ) -> EmpiricalDistribution:
        """Generate typical properties of the given type."""
        typical_property_sets = [
            infinite_network.calc_base_property_value_set(property_type)
            for infinite_network in infinite_networks
        ]

        distribution = EmpiricalDistribution([
            value for property_set in typical_property_sets for value in property_set
        ])
        return distribution

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
