#!/usr/bin/env python3
"""Factory methods for creating theoretical distributions."""

from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.normal_distribution import NormalDistribution
from distribution.theoretical.poisson_distribution import PoissonDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from distribution.theoretical.stable_distribution import StableDistribution
from distribution.theoretical.uniform_distribution import UniformDistribution


def create_theoretical_distribution(distribution_type: TheoreticalDistribution.Type) -> TheoreticalDistribution:
    """Create a given theoretical distributions."""
    if distribution_type == TheoreticalDistribution.Type.NORMAL:
        return NormalDistribution()
    if distribution_type == TheoreticalDistribution.Type.POISSON:
        return PoissonDistribution()
    if distribution_type == TheoreticalDistribution.Type.POWER_LAW:
        return PowerLawDistribution()
    if distribution_type == TheoreticalDistribution.Type.STABLE:
        return StableDistribution()
    if distribution_type == TheoreticalDistribution.Type.UNIFORM:
        return UniformDistribution()
    raise NotImplementedError


def create_distribution_parameters(
    distribution_type: TheoreticalDistribution.Type
) -> TheoreticalDistribution.Parameters:
    """Create a given theoretical distributions."""
    if distribution_type == TheoreticalDistribution.Type.NORMAL:
        return NormalDistribution.Parameters()
    if distribution_type == TheoreticalDistribution.Type.POISSON:
        return PoissonDistribution.Parameters()
    if distribution_type == TheoreticalDistribution.Type.POWER_LAW:
        return PowerLawDistribution.Parameters()
    if distribution_type == TheoreticalDistribution.Type.STABLE:
        return StableDistribution.Parameters()
    if distribution_type == TheoreticalDistribution.Type.UNIFORM:
        return UniformDistribution.Parameters()
    raise NotImplementedError


def create_fitting_method(
    distribution_type: TheoreticalDistribution.Type
) -> TheoreticalDistribution.FittingMethod:
    """Create a given theoretical distributions."""
    if distribution_type == TheoreticalDistribution.Type.NORMAL:
        return NormalDistribution.FittingMethod.MAXIMUM_LIKELIHOOD
    if distribution_type == TheoreticalDistribution.Type.POISSON:
        return PoissonDistribution.FittingMethod.MAXIMUM_LIKELIHOOD
    if distribution_type == TheoreticalDistribution.Type.POWER_LAW:
        return PowerLawDistribution.FittingMethod.MAXIMUM_LIKELIHOOD_DETERMINISTIC_DOMAIN
    if distribution_type == TheoreticalDistribution.Type.STABLE:
        return StableDistribution.FittingMethod.MLE_SCIPY
    if distribution_type == TheoreticalDistribution.Type.UNIFORM:
        return UniformDistribution.FittingMethod.DEFAULT
    raise NotImplementedError


def create_default_fitting_parameters(
    distribution_type: TheoreticalDistribution.Type
) -> TheoreticalDistribution.FittingParameters:
    """Create a given theoretical distributions."""
    default_parameters = create_distribution_parameters(distribution_type)
    fitting_method = create_fitting_method(distribution_type)
    if distribution_type == TheoreticalDistribution.Type.NORMAL:
        return NormalDistribution.FittingParameters(default_parameters, fitting_method)
    if distribution_type == TheoreticalDistribution.Type.POISSON:
        return PoissonDistribution.FittingParameters(default_parameters, fitting_method)
    if distribution_type == TheoreticalDistribution.Type.POWER_LAW:
        return PowerLawDistribution.FittingParameters(default_parameters, fitting_method)
    if distribution_type == TheoreticalDistribution.Type.STABLE:
        return StableDistribution.FittingParameters(default_parameters, fitting_method)
    if distribution_type == TheoreticalDistribution.Type.UNIFORM:
        return UniformDistribution.FittingParameters(default_parameters, fitting_method)
    raise NotImplementedError
