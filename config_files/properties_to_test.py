#!/usr/bin/env python3
"""Specification of scalar property params to test."""

import numpy as np

from config_files.distribution_fitting_params import (
    POWER_LAW_FITTING_MINIMUM_VALUE_DATA,
    POWER_LAW_FITTING_MINIMUM_VALUE_MODEL,
)
from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import (
    create_fitting_parameters_normal,
    create_power_law_fitting_parameters,
)
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.normal_distribution import NormalDistribution
from distribution.theoretical.poisson_distribution import PoissonDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from distribution.theoretical.stable_distribution import StableDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty


FINITE_SCALAR_PROPERTY_NAMES_TO_TEST = [
    'vertex_degree_exponent',
    'edge_degree_exponent',
    'vertex_interaction_degree_exponent',
    'interaction_vertex_degree_exponent',
    # 'triangle_degree_exponent',
    # 'average_interaction_degree_normal_mle',
    # 'num_of_isolated_vertices_normal_mle',
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


INFINITE_SCALAR_PROPERTY_NAMES_TO_TEST = [
    'vertex_degree_exponent',
    'edge_degree_exponent',
    'vertex_interaction_degree_exponent',
    'interaction_vertex_degree_exponent',
    # 'triangle_degree_exponent',
    'average_interaction_degree_normal_mle',
    # 'num_of_edges_normal_mle',
    # 'num_of_edges_normal_match_quantile',
    # 'num_of_edges_stable',
    # 'num_of_triangles_normal_mle',
    # 'num_of_triangles_normal_match_quantile',
    # 'num_of_triangles_stable',
]


def _get_power_law_exponent_model(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_params = create_power_law_fitting_parameters(POWER_LAW_FITTING_MINIMUM_VALUE_MODEL)
    approximation.fit(fitting_params)
    # print(approximation.info()['empirical_values_in_theoretical_domain'])
    assert isinstance(approximation.theoretical, PowerLawDistribution)
    return approximation.theoretical.parameters.exponent


def _get_power_law_exponent_data_set(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(
        empirical_distribution,
        TheoreticalDistribution.Type.POWER_LAW
    )
    fitting_params = create_power_law_fitting_parameters(POWER_LAW_FITTING_MINIMUM_VALUE_DATA)
    approximation.fit(fitting_params)
    assert isinstance(approximation.theoretical, PowerLawDistribution)
    return approximation.theoretical.parameters.exponent


def _get_power_law_exponent_infinite_model_model(list_of_degree_sets: list[list[int]]) -> float:
    flat_list = [degree for degree_set in list_of_degree_sets for degree in degree_set]
    return _get_power_law_exponent_model(EmpiricalDistribution(flat_list))


def _get_power_law_exponent_infinite_model_data_set(list_of_degree_sets: list[list[int]]) -> float:
    flat_list = [degree for degree_set in list_of_degree_sets for degree in degree_set]
    return _get_power_law_exponent_data_set(EmpiricalDistribution(flat_list))


def _get_poisson_parameter(empirical_distribution: EmpiricalDistribution) -> float:
    approximation = DistributionApproximation(empirical_distribution, TheoreticalDistribution.Type.POISSON)
    approximation.fit(create_fitting_parameters_normal())
    assert isinstance(approximation.theoretical, PoissonDistribution)
    return approximation.theoretical.parameters.lambda_


def _get_mean(empirical_distribution: EmpiricalDistribution) -> float:
    return empirical_distribution.value_sequence.mean()


def _get_mean_infinite(list_of_value_sets: list[list[int]]) -> float:
    return np.nanmean([value for set_ in list_of_value_sets for value in set_])


def get_finite_scalar_property_params(gamma: float) -> list[DerivedNetworkProperty]:
    """Return the scalar property params to test for finite networks."""
    all_properties = (
        DerivedNetworkProperty(name='vertex_degree_exponent',
                               source_base_property=BaseNetworkProperty.vertex_edge_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_model,
                               calculator_data_set=_get_power_law_exponent_data_set),
        DerivedNetworkProperty(name='edge_degree_exponent',
                               source_base_property=BaseNetworkProperty.edge_triangle_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_model,
                               calculator_data_set=_get_power_law_exponent_data_set),
        DerivedNetworkProperty(name='triangle_degree_exponent',
                               source_base_property=BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_model,
                               calculator_data_set=_get_power_law_exponent_data_set),
        DerivedNetworkProperty(name='average_interaction_degree_normal_mle',
                               source_base_property=BaseNetworkProperty.vertex_interaction_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_mean,
                               calculator_data_set=_get_mean),
        DerivedNetworkProperty(name='vertex_interaction_degree_exponent',
                               source_base_property=BaseNetworkProperty.vertex_interaction_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_model,
                               calculator_data_set=_get_power_law_exponent_data_set),
        DerivedNetworkProperty(name='interaction_vertex_degree_exponent',
                               source_base_property=BaseNetworkProperty(
                                   BaseNetworkProperty.interaction_vertex_degree_distribution
                                ),
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_model,
                               calculator_data_set=_get_power_law_exponent_data_set),
        DerivedNetworkProperty(name='num_of_isolated_vertices_normal_mle',
                               source_base_property=BaseNetworkProperty(
                                   BaseNetworkProperty.vertex_interaction_degree_distribution),
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=lambda degree_distribution:
                                   degree_distribution.calc_value_counts()[0, 1]),
        DerivedNetworkProperty(name='num_of_edges_normal_mle',
                               source_base_property=BaseNetworkProperty.num_of_edges,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal()),
        DerivedNetworkProperty(name='num_of_edges_normal_match_quantile',
                               source_base_property=BaseNetworkProperty.num_of_edges,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=NormalDistribution.FittingParameters(
                                   NormalDistribution.DomainCalculation(),
                                   NormalDistribution.ParameterFittingMatchQuantile(
                                       NormalDistribution.ParameterFitting.Method.MATCH_QUANTILE,
                                       NormalDistribution.DistributionParameters(),
                                       quantile=0.025,
                                   ),
                               )),
        DerivedNetworkProperty(name='num_of_edges_stable',
                               source_base_property=BaseNetworkProperty.num_of_edges,
                               theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
                               fitting_parameters=StableDistribution.FittingParameters(
                                   StableDistribution.DomainCalculation(),
                                   StableDistribution.ParameterFitting(
                                       StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                                       StableDistribution.DistributionParameters(
                                           alpha=min(1 / gamma, 2.),
                                           # alpha=np.nan,
                                           beta=1.,
                                           location=np.nan,
                                           scale=np.nan,
                                       ),
                                   ),
                               )),
        DerivedNetworkProperty(name='num_of_triangles_normal_mle',
                               source_base_property=BaseNetworkProperty.num_of_triangles,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal()),
        DerivedNetworkProperty(name='num_of_triangles_normal_match_quantile',
                               source_base_property=BaseNetworkProperty.num_of_triangles,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=NormalDistribution.FittingParameters(
                                   NormalDistribution.DomainCalculation(),
                                   NormalDistribution.ParameterFittingMatchQuantile(
                                       NormalDistribution.ParameterFitting.Method.MATCH_QUANTILE,
                                       NormalDistribution.DistributionParameters(),
                                       quantile=0.025,
                                   ),
                               )),
        DerivedNetworkProperty(name='num_of_triangles_stable',
                               source_base_property=BaseNetworkProperty.num_of_triangles,
                               theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
                               fitting_parameters=StableDistribution.FittingParameters(
                                   StableDistribution.DomainCalculation(),
                                   StableDistribution.ParameterFitting(
                                       StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                                       StableDistribution.DistributionParameters(
                                           alpha=min(1 / gamma, 2.),
                                           # alpha=np.nan,
                                           beta=1.,
                                           location=np.nan,
                                           scale=np.nan
                                       ),
                                   ),
                               )),
        DerivedNetworkProperty(name='betti_number_0_normal',
                               source_base_property=BaseNetworkProperty.betti_numbers,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=lambda betti_numbers: betti_numbers[0]),
        DerivedNetworkProperty(name='betti_number_0_stable',
                               source_base_property=BaseNetworkProperty.betti_numbers,
                               theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
                               fitting_parameters=StableDistribution.FittingParameters(
                                   StableDistribution.DomainCalculation(),
                                   StableDistribution.ParameterFitting(
                                       StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                                       StableDistribution.DistributionParameters(
                                           alpha=min(1 / gamma, 2.),
                                           # alpha=np.nan,
                                           beta=-1.,
                                           location=np.nan,
                                           scale=np.nan
                                       ),
                                   ),
                               ),
                               calculator_default=lambda betti_numbers: betti_numbers[0]),
        DerivedNetworkProperty(name='betti_number_1_normal',
                               source_base_property=BaseNetworkProperty.betti_numbers,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=lambda betti_numbers: betti_numbers[1]),
        DerivedNetworkProperty(name='betti_number_1_stable',
                               source_base_property=BaseNetworkProperty.betti_numbers,
                               theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
                               fitting_parameters=StableDistribution.FittingParameters(
                                   StableDistribution.DomainCalculation(),
                                   StableDistribution.ParameterFitting(
                                       StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                                       StableDistribution.DistributionParameters(
                                           alpha=min(1 / gamma, 2.),
                                           # alpha=np.nan,
                                           beta=-1.,
                                           location=np.nan,
                                           scale=np.nan
                                       ),
                                   ),
                               ),
                               calculator_default=lambda betti_numbers: betti_numbers[1]),
        DerivedNetworkProperty(name='betti_number_2_normal',
                               source_base_property=BaseNetworkProperty.betti_numbers,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=lambda betti_numbers: betti_numbers[2]),
        DerivedNetworkProperty(name='betti_number_2_stable',
                               source_base_property=BaseNetworkProperty.betti_numbers,
                               theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
                               fitting_parameters=StableDistribution.FittingParameters(
                                   StableDistribution.DomainCalculation(),
                                   StableDistribution.ParameterFitting(
                                       StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                                       StableDistribution.DistributionParameters(
                                           alpha=min(1 / gamma, 2.),
                                           beta=-1.,
                                           location=np.nan,
                                           scale=np.nan
                                       ),
                                   ),
                               ),
                               calculator_default=lambda betti_numbers: betti_numbers[2]),
    )
    return [property_ for property_ in all_properties if property_.name in FINITE_SCALAR_PROPERTY_NAMES_TO_TEST]


def get_infinite_scalar_property_params(gamma: float) -> list[DerivedNetworkProperty]:
    """Return the scalar property params to test for infinite networks."""
    all_properties = (
        DerivedNetworkProperty(name='vertex_degree_exponent',
                               source_base_property=BaseNetworkProperty.vertex_edge_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_infinite_model_model,
                               calculator_data_set=_get_power_law_exponent_infinite_model_data_set),
        DerivedNetworkProperty(name='edge_degree_exponent',
                               source_base_property=BaseNetworkProperty.edge_triangle_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_infinite_model_model,
                               calculator_data_set=_get_power_law_exponent_infinite_model_data_set),
        DerivedNetworkProperty(name='triangle_degree_exponent',
                               source_base_property=BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_infinite_model_model,
                               calculator_data_set=_get_power_law_exponent_infinite_model_data_set),
        DerivedNetworkProperty(name='average_interaction_degree_normal_mle',
                               source_base_property=BaseNetworkProperty.vertex_interaction_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_mean_infinite,
                               calculator_data_set=_get_mean_infinite),
        DerivedNetworkProperty(name='vertex_interaction_degree_exponent',
                               source_base_property=BaseNetworkProperty.vertex_interaction_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_infinite_model_model,
                               calculator_data_set=_get_power_law_exponent_infinite_model_data_set),
        DerivedNetworkProperty(name='interaction_vertex_degree_exponent',
                               source_base_property=BaseNetworkProperty.interaction_vertex_degree_distribution,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=_get_power_law_exponent_infinite_model_model,
                               calculator_data_set=_get_power_law_exponent_infinite_model_data_set),
        DerivedNetworkProperty(name='num_of_edges_normal_mle',
                               source_base_property=BaseNetworkProperty.num_of_edges,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=sum,),
        DerivedNetworkProperty(name='num_of_edges_normal_match_quantile',
                               source_base_property=BaseNetworkProperty.num_of_edges,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=NormalDistribution.FittingParameters(
                                   NormalDistribution.DomainCalculation(),
                                   NormalDistribution.ParameterFittingMatchQuantile(
                                       NormalDistribution.ParameterFitting.Method.MATCH_QUANTILE,
                                       NormalDistribution.DistributionParameters(),
                                       quantile=0.025,
                                   ),
                               ),
                               calculator_default=sum,),
        DerivedNetworkProperty(name='num_of_edges_stable',
                               source_base_property=BaseNetworkProperty.num_of_edges,
                               theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
                               fitting_parameters=StableDistribution.FittingParameters(
                                   StableDistribution.DomainCalculation(),
                                   StableDistribution.ParameterFitting(
                                       StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                                       StableDistribution.DistributionParameters(
                                           alpha=min(1 / gamma, 2.),
                                           # alpha=np.nan,
                                           beta=1.,
                                           location=np.nan,
                                           scale=np.nan
                                       ),
                                   ),
                               ),
                               calculator_default=sum,),
        DerivedNetworkProperty(name='num_of_triangles_normal_mle',
                               source_base_property=BaseNetworkProperty.num_of_triangles,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=create_fitting_parameters_normal(),
                               calculator_default=sum),
        DerivedNetworkProperty(name='num_of_triangles_normal_match_quantile',
                               source_base_property=BaseNetworkProperty.num_of_triangles,
                               theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
                               fitting_parameters=NormalDistribution.FittingParameters(
                                   NormalDistribution.DomainCalculation(),
                                   NormalDistribution.ParameterFittingMatchQuantile(
                                       NormalDistribution.ParameterFitting.Method.MATCH_QUANTILE,
                                       NormalDistribution.DistributionParameters(),
                                       quantile=0.025,
                                   ),
                               ),
                               calculator_default=sum),
        DerivedNetworkProperty(name='num_of_triangles_stable',
                               source_base_property=BaseNetworkProperty.num_of_triangles,
                               theoretical_approximation_type=TheoreticalDistribution.Type.STABLE,
                               fitting_parameters=StableDistribution.FittingParameters(
                                   StableDistribution.DomainCalculation(),
                                   StableDistribution.ParameterFitting(
                                       StableDistribution.ParameterFitting.Method.MLE_SCIPY,
                                       StableDistribution.DistributionParameters(
                                           alpha=min(1 / gamma, 2.),
                                           # alpha=np.nan,
                                           beta=1.,
                                           location=np.nan,
                                           scale=np.nan
                                       ),
                                   ),
                               ),
                               calculator_default=sum),
    )
    return [property_ for property_ in all_properties if property_.name in INFINITE_SCALAR_PROPERTY_NAMES_TO_TEST]
