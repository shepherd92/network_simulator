#!/usr/bin/env python3
"""Specification of scalar property params to test."""

import numpy as np

from config_files.model_config.age_dependent_random_simplex_config import AGE_DEPENDENT_RANDOM_SIMPLEX_MODEL_PARAMETERS
from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.normal_distribution import NormalDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from distribution.theoretical.stable_distribution import StableDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty


def _get_power_law_exponent(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(empirical_distribution, TheoreticalDistribution.Type.POWER_LAW)
    approximation.fit(PowerLawDistribution.FittingParameters(
        PowerLawDistribution.Parameters(),
        PowerLawDistribution.FittingMethod.MAXIMUM_LIKELIHOOD_MLE_DOMAIN,
    ))
    assert isinstance(approximation.theoretical, PowerLawDistribution)
    return approximation.theoretical.parameters.exponent


GAMMA = AGE_DEPENDENT_RANDOM_SIMPLEX_MODEL_PARAMETERS.gamma


SCALAR_PROPERTY_PARAMS_TO_TEST: tuple[DerivedNetworkProperty, ...] = (
    DerivedNetworkProperty(
        name='num_of_edges_normal_mle',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            fixed_parameters=NormalDistribution.Parameters(),
            fitting_method=NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
    ),
    DerivedNetworkProperty(
        name='num_of_edges_normal_match_quantile',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            fixed_parameters=NormalDistribution.Parameters(),
            fitting_method=NormalDistribution.FittingMethod.MATCH_2_DOT_5_PERCENTILE,
        ),
    ),
    DerivedNetworkProperty(
        name='num_of_edges_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.Parameters(alpha=1 / GAMMA, beta=1., location=np.nan, scale=np.nan),
            StableDistribution.FittingMethod.MLE_SCIPY,
        ),
    ),
    DerivedNetworkProperty(
        name='vertex_degree_exponent_finite',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION,
            BaseNetworkProperty.CalculationMethod.NETWORK,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.Parameters(),
            NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
        calculator=_get_power_law_exponent
    ),
    DerivedNetworkProperty(
        name='vertex_degree_exponent_infinite',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION,
            BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.Parameters(),
            NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
        calculator=_get_power_law_exponent
    ),
    DerivedNetworkProperty(
        name='edge_degree_exponent_finite',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1,
            BaseNetworkProperty.CalculationMethod.NETWORK,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.Parameters(),
            NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
        calculator=lambda distribution: _get_power_law_exponent(distribution)
    ),
    DerivedNetworkProperty(
        name='edge_degree_exponent_infinite',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1,
            BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.Parameters(),
            NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
        calculator=lambda distribution: _get_power_law_exponent(distribution)
    ),
    DerivedNetworkProperty(
        name='betti_number_0_normal',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.Parameters(),
            NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
        calculator=lambda betti_numbers: betti_numbers[0, 1],
    ),
    DerivedNetworkProperty(
        name='betti_number_0_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.Parameters(),
            StableDistribution.FittingMethod.MLE_SCIPY,
        ),
        calculator=lambda betti_numbers: betti_numbers[0, 1],
    ),
    DerivedNetworkProperty(
        name='betti_number_1_normal',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.Parameters(),
            NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
        calculator=lambda betti_numbers: betti_numbers[1, 1],
    ),
    DerivedNetworkProperty(
        name='betti_number_1_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.Parameters(),
            StableDistribution.FittingMethod.MLE_SCIPY,
        ),
        calculator=lambda betti_numbers: betti_numbers[1, 1],
    ),
)
