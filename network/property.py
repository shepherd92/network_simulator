#!/usr/bin/env python3
"""This module is responsible for calculating the summaries for a given network."""

from dataclasses import dataclass
from enum import Enum, auto, unique
from typing import Any, Callable

from distribution.approximation import DistributionApproximation
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


@unique
class BaseNetworkProperty(Enum):
    """Represent the calculatable properties of a model."""
    num_of_vertices: int = auto()
    num_of_edges: int = auto()
    num_of_triangles: int = auto()
    num_of_interactions: int = auto()
    edges: int = auto()
    mean_degree: int = auto()
    max_degree: int = auto()
    interaction_vertex_degree_distribution: int = auto()
    simplex_dimension_distribution: int = auto()
    vertex_edge_degree_distribution: int = auto()
    vertex_interaction_degree_distribution: int = auto()
    edge_interaction_degree_distribution: int = auto()
    triangle_interaction_degree_distribution: int = auto()
    edge_triangle_degree_distribution: int = auto()
    triangle_tetrahedra_degree_distribution: int = auto()
    betti_numbers: int = auto()
    betti_numbers_by_component: int = auto()
    num_of_vertices_by_component: int = auto()
    persistence_intervals: int = auto()
    persistence_pairs: int = auto()


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
    calculator_dataset: Callable[[Any], float | int] | None = None
    directly_calculated_from_model: bool = False


@dataclass
class ScalarNetworkPropertyReport:
    """Represent everything about a scalar network property for reporting."""

    params: DerivedNetworkProperty
    distributions: DistributionApproximation
    test_results: DistributionApproximation.TestResult
    data_point: float
