#!/usr/bin/env python3
"""This module is responsible for calculating the summaries for a given network."""

from dataclasses import dataclass
from enum import Enum, auto
from typing import Any, Callable

from distribution.approximation import DistributionApproximation
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


@dataclass
class BaseNetworkProperty:
    """Represent the calculatable properties of a model."""

    class Type(Enum):
        """Represent the calculatable property types of a network."""

        NUM_OF_NODES: int = auto()
        NUM_OF_EDGES: int = auto()
        NUM_OF_TRIANGLES: int = auto()
        AVERAGE_DEGREE: int = auto()
        MAX_DEGREE: int = auto()
        AVG_CLUSTERING: int = auto()
        NUM_OF_CONNECTED_COMPONENTS: int = auto()
        DIMENSION: int = auto()
        NUM_OF_SIMPLICES: int = auto()
        INTERACTION_DIMENSION_DISTRIBUTION: int = auto()
        SIMPLEX_DIMENSION_DISTRIBUTION: int = auto()
        FACET_DIMENSION_DISTRIBUTION: int = auto()
        DEGREE_DISTRIBUTION: int = auto()
        IN_DEGREE_DISTRIBUTION: int = auto()
        OUT_DEGREE_DISTRIBUTION: int = auto()
        HIGHER_ORDER_DEGREE_DISTRIBUTION_1: int = auto()
        HIGHER_ORDER_DEGREE_DISTRIBUTION_2: int = auto()
        HIGHER_ORDER_DEGREE_DISTRIBUTION_3: int = auto()
        BETTI_NUMBERS: int = auto()
        BETTI_NUMBERS_BY_COMPONENT: int = auto()
        VERTICES_BY_COMPONENT: int = auto()
        PERSISTENCE: int = auto()

    class CalculationMethod(Enum):
        """Represent the calculatable property types of a network."""

        NETWORK: int = auto()
        TYPICAL_OBJECT: int = auto()

    property_type: Type
    calculation_method: CalculationMethod = CalculationMethod.NETWORK


@dataclass
class DerivedNetworkProperty:
    """Represent a scalar derived network property that can be used for fitting and testing.

    The calculator member receives the calculated possible base network property,
    and extracts the scalar value of interest.
    """

    name: str
    source_base_property: BaseNetworkProperty
    theoretical_approximation_type: TheoreticalDistribution.Type
    fitting_parameters: TheoreticalDistribution.FittingParameters
    calculator_default: Callable[[Any], float | int] = lambda property_value: property_value
    calculator_data_set: Callable[[Any], float | int] | None = None


@dataclass
class ScalarNetworkPropertyReport:
    """Represent everything about a scalar network property for reporting."""

    params: DerivedNetworkProperty
    distributions: DistributionApproximation
    test_results: DistributionApproximation.TestResult
    data_point: float
