#!/usr/bin/env python3
"""Specification of scalar property params to test."""

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty


def _get_power_law_exponent(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(empirical_distribution, TheoreticalDistribution.Type.POWER_LAW)
    approximation.fit()
    assert isinstance(approximation.theoretical, PowerLawDistribution)
    return approximation.theoretical.parameters.exponent


SCALAR_PROPERTY_PARAMS_TO_TEST: tuple[DerivedNetworkProperty, ...] = (
    DerivedNetworkProperty(
        name='num_of_edges',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
    ),
    DerivedNetworkProperty(
        name='num_of_edges_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
    ),
    DerivedNetworkProperty(
        name='vertex_degree_exponent',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION,
            BaseNetworkProperty.CalculationMethod.NETWORK,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        calculator=_get_power_law_exponent
    ),
    DerivedNetworkProperty(
        name='edge_degree_exponent',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTIONS,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        calculator=lambda distributions: _get_power_law_exponent(distributions[1])
    ),
    DerivedNetworkProperty(
        name='triangle_degree_exponent',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTIONS,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        calculator=lambda distributions: _get_power_law_exponent(distributions[2])
    ),
    DerivedNetworkProperty(
        name='betti_number_0',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        calculator=lambda betti_numbers: betti_numbers[0],
    ),
    DerivedNetworkProperty(
        name='betti_number_1',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        calculator=lambda betti_numbers: betti_numbers[1],
    ),
    DerivedNetworkProperty(
        name='betti_number_2',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        calculator=lambda betti_numbers: betti_numbers[2],
    ),
)
