#!/usr/bin/env python3
"""Specification of scalar property params to test."""

from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty


def _get_power_law_exponent(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(empirical_distribution, TheoreticalDistribution.Type.POWER_LAW)
    assert isinstance(approximation.theoretical, PowerLawDistribution)
    return approximation.theoretical.parameters.exponent


SCALAR_PROPERTY_PARAMS_TO_TEST: tuple[DerivedNetworkProperty, ...] = (
    DerivedNetworkProperty(
        name='num_of_nodes',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_NODES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
    ),
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
        name='in_degree_exponent',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION,
            BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        calculator=_get_power_law_exponent
    ),
    DerivedNetworkProperty(
        name='betti_number_1',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        calculator=lambda betti_numbers: betti_numbers[1],
    ),
    # ScalarNetworkPropertyParams(
    #     name='node_degree_distribution_power_law_exponent',
    #     source_base_property_type=BaseNetworkPropertyType.ORDINARY_DEGREE_DISTRIBUTIONS,
    #     calculator=lambda property_value: property_value['in'].theoretical.exponent,
    #     theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
    # ),
    # ScalarNetworkPropertyParams(
    #     name='edge_degree_distribution_power_law_exponent',
    #     source_base_property_type=BaseNetworkPropertyType.HIGHER_ORDER_DEGREE_DISTRIBUTIONS,
    #     calculator=lambda property_value: property_value[1].theoretical.exponent,
    #     theoretical_approximation=NormalDistribution(),
    # ),
    # ScalarNetworkPropertyParams(
    #     name='betti_number_1',
    #     source_base_property_type=BaseNetworkPropertyType.BETTI_NUMBERS,
    #     calculator=lambda property_value: property_value[1],
    #     theoretical_approximation=NormalDistribution(),
    # ),
    # ScalarNetworkPropertyParams(
    #     name='degree_distribution_power_law_hypothesis_ks_statistic',
    #     source_base_property_type=BaseNetworkPropertyType.POWER_LAW_PARAMETERS,
    #     calculator=lambda property_value: property_value.kolmogorov_smirnov_statistics,
    #     theoretical_approximation=NormalDistribution(),
    # ),
    # ScalarNetworkPropertyParams(
    #     name='avg_clustering',
    #     source_base_property_type=BaseNetworkPropertyType.AVG_CLUSTERING,
    #     theoretical_approximation=NormalDistribution(),
    # ),
    # ScalarNetworkPropertyParams(
    #     name='facet_dimension_distribution_dim_1',
    #     source_base_property_type=BaseNetworkPropertyType.FACET_DIMENSION_DISTRIBUTION,  # noqa: E251
    #     calculator=lambda property_value: property_value[1],
    #     theoretical_approximation=NormalDistribution(),
    # ),
)
