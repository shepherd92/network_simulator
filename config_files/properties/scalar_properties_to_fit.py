#!/usr/bin/env python3
"""Specification of scalar property params to fit."""

from distribution.theoretical.normal_distribution import NormalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty


SCALAR_PROPERTY_PARAMS_TO_FIT: tuple[DerivedNetworkProperty, ...] = (
    DerivedNetworkProperty(
        name='num_of_edges',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        fitting_parameters=NormalDistribution.FittingParameters(
            fixed_parameters=NormalDistribution.Parameters(),
            fitting_method=NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD,
        ),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
    ),
)
