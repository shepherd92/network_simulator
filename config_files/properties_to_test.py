#!/usr/bin/env python3
"""Specification of scalar property params to test."""

import numpy as np

from config_files.model_config import HYPERGRAPH_MODEL_PARAMETERS
from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import (
    create_fitting_parameters_normal,
    create_fitting_parameters_power_law_model,
    create_fitting_parameters_power_law_data_set,
)
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.normal_distribution import NormalDistribution
from distribution.theoretical.poisson_distribution import PoissonDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from distribution.theoretical.stable_distribution import StableDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty


SCALAR_PROPERTY_NAMES_TO_TEST = [
    'vertex_degree_exponent',
    'edge_degree_exponent',
    # 'triangle_degree_exponent',
    # 'average_interaction_degree_normal_mle',
    # 'num_of_edges_normal_mle',
    # 'num_of_edges_normal_match_quantile',
    # 'num_of_edges_stable',
    # 'num_of_triangles_normal_mle',
    # 'num_of_triangles_normal_match_quantile',
    # 'num_of_triangles_stable',
    # 'betti_number_0_normal',
    # 'betti_number_0_stable',
    # 'betti_number_1_normal',
    # 'betti_number_1_stable',
    # 'betti_number_2_normal',
    # 'betti_number_2_stable',
]

GAMMA = HYPERGRAPH_MODEL_PARAMETERS.gamma


def _get_power_law_exponent_model(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_params = create_fitting_parameters_power_law_model()
    approximation.fit(fitting_params)
    assert isinstance(approximation.theoretical, PowerLawDistribution)
    return approximation.theoretical.parameters.exponent


def _get_power_law_exponent_data_set(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_params = create_fitting_parameters_power_law_data_set()
    approximation.fit(fitting_params)
    assert isinstance(approximation.theoretical, PowerLawDistribution)
    return approximation.theoretical.parameters.exponent


def _get_poisson_parameter(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(empirical_distribution, TheoreticalDistribution.Type.POISSON)
    approximation.fit(create_fitting_parameters_normal())
    assert isinstance(approximation.theoretical, PoissonDistribution)
    return approximation.theoretical.parameters.lambda_


def _get_mean(empirical_distribution: EmpiricalDistribution) -> float:
    return empirical_distribution.value_sequence.mean()


ALL_SCALAR_PROPERTY_PARAMS: tuple[DerivedNetworkProperty, ...] = (

    DerivedNetworkProperty(
        name='vertex_degree_exponent',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION,
            BaseNetworkProperty.CalculationMethod.NETWORK,
            # BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
        calculator_default=_get_power_law_exponent_model,
        calculator_data_set=_get_power_law_exponent_data_set,
    ),
    DerivedNetworkProperty(
        name='edge_degree_exponent',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1,
            BaseNetworkProperty.CalculationMethod.NETWORK,
            # BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
        calculator_default=_get_power_law_exponent_model,
        calculator_data_set=_get_power_law_exponent_data_set,
    ),
    DerivedNetworkProperty(
        name='triangle_degree_exponent',
        source_base_property=BaseNetworkProperty(
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2,
            # BaseNetworkProperty.CalculationMethod.NETWORK,
            BaseNetworkProperty.CalculationMethod.TYPICAL_OBJECT,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
        calculator_default=_get_power_law_exponent_model,
        calculator_data_set=_get_power_law_exponent_data_set,
    ),
    DerivedNetworkProperty(
        name='average_interaction_degree_normal_mle',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.VERTEX_INTERACTION_DEGREE_DISTRIBUTION),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
        calculator_default=_get_mean,
        calculator_data_set=_get_mean,
    ),
    DerivedNetworkProperty(
        name='num_of_edges_normal_mle',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
    ),
    DerivedNetworkProperty(
        name='num_of_edges_normal_match_quantile',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.DomainCalculation(),
            NormalDistribution.ParameterFittingMatchQuantile(
                NormalDistribution.ParameterFitting.Method.MATCH_QUANTILE,
                NormalDistribution.Parameters(),
                quantile=0.025,
            ),
        ),
    ),
    DerivedNetworkProperty(
        name='num_of_edges_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.DomainCalculation(),
            StableDistribution.ParameterFitting(
                StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                StableDistribution.Parameters(
                    alpha=min(1 / GAMMA, 2.),
                    # alpha=np.nan,
                    beta=1.,
                    location=np.nan,
                    scale=np.nan
                ),
            ),
        ),
    ),
    DerivedNetworkProperty(
        name='num_of_triangles_normal_mle',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_TRIANGLES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
    ),
    DerivedNetworkProperty(
        name='num_of_triangles_normal_match_quantile',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_TRIANGLES),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=NormalDistribution.FittingParameters(
            NormalDistribution.DomainCalculation(),
            NormalDistribution.ParameterFittingMatchQuantile(
                NormalDistribution.ParameterFitting.Method.MATCH_QUANTILE,
                NormalDistribution.Parameters(),
                quantile=0.025,
            ),
        ),
    ),
    DerivedNetworkProperty(
        name='num_of_triangles_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_TRIANGLES),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.DomainCalculation(),
            StableDistribution.ParameterFitting(
                StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                StableDistribution.Parameters(
                    alpha=min(1 / GAMMA, 2.),
                    # alpha=np.nan,
                    beta=1.,
                    location=np.nan,
                    scale=np.nan
                ),
            ),
        ),
    ),
    DerivedNetworkProperty(
        name='betti_number_0_normal',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
        calculator_default=lambda betti_numbers: betti_numbers[0],
    ),
    DerivedNetworkProperty(
        name='betti_number_0_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.DomainCalculation(),
            StableDistribution.ParameterFitting(
                StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                StableDistribution.Parameters(
                    alpha=min(1 / GAMMA, 2.),
                    # alpha=np.nan,
                    beta=-1.,
                    location=np.nan,
                    scale=np.nan
                ),
            ),
        ),
        calculator_default=lambda betti_numbers: betti_numbers[0],
    ),
    DerivedNetworkProperty(
        name='betti_number_1_normal',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
        calculator_default=lambda betti_numbers: betti_numbers[1],
    ),
    DerivedNetworkProperty(
        name='betti_number_1_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.DomainCalculation(),
            StableDistribution.ParameterFitting(
                StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                StableDistribution.Parameters(
                    alpha=min(1 / GAMMA, 2.),
                    # alpha=np.nan,
                    beta=-1.,
                    location=np.nan,
                    scale=np.nan
                ),
            ),
        ),
        calculator_default=lambda betti_numbers: betti_numbers[1],
    ),
    DerivedNetworkProperty(
        name='betti_number_2_normal',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
        fitting_parameters=create_fitting_parameters_normal(),
        calculator_default=lambda betti_numbers: betti_numbers[2],
    ),
    DerivedNetworkProperty(
        name='betti_number_2_stable',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.BETTI_NUMBERS),
        theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
        fitting_parameters=StableDistribution.FittingParameters(
            StableDistribution.DomainCalculation(),
            StableDistribution.ParameterFitting(
                StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                StableDistribution.Parameters(
                    alpha=min(1 / GAMMA, 2.),
                    beta=-1.,
                    location=np.nan,
                    scale=np.nan
                ),
            ),
        ),
        calculator_default=lambda betti_numbers: betti_numbers[2],
    ),
)


SCALAR_PROPERTY_PARAMS_TO_TEST: list[DerivedNetworkProperty] = [
    property_
    for property_ in ALL_SCALAR_PROPERTY_PARAMS
    if property_.name in SCALAR_PROPERTY_NAMES_TO_TEST
]
