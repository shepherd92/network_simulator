#!/usr/bin/env python3
"""Factory methods for creating theoretical distributions."""

import numpy as np

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


def create_fitting_parameters_normal() -> NormalDistribution.FittingParameters:
    """Create fitting parameters for normal distribution."""
    return NormalDistribution.FittingParameters(
        NormalDistribution.DomainCalculation(),
        NormalDistribution.ParameterFitting(
            NormalDistribution.ParameterFitting.Method.MAXIMUM_LIKELIHOOD,
            NormalDistribution.DistributionParameters(),
        )
    )


def create_fitting_parameters_poisson() -> PoissonDistribution.FittingParameters:
    """Create fitting parameters for normal distribution."""
    return PoissonDistribution.FittingParameters(
        PoissonDistribution.DomainCalculation(),
        PoissonDistribution.ParameterFitting(
            PoissonDistribution.ParameterFitting.Method.MAXIMUM_LIKELIHOOD,
            PoissonDistribution.DistributionParameters(),
        )
    )


def create_power_law_fitting_parameters(minimum_value: float) -> PowerLawDistribution.FittingParameters:
    """Create fitting parameters for normal distribution."""
    return PowerLawDistribution.FittingParameters(
        PowerLawDistribution.DeterministicDomain(
            PowerLawDistribution.DomainCalculation.Method.DETERMINISTIC,
            min_=minimum_value,
            max_=np.inf,
        ),
        # PowerLawDistribution.MaximumLikelihoodDomain(
        #     PowerLawDistribution.DomainCalculation.Method.MAXIMUM_LIKELIHOOD,
        # ),
        PowerLawDistribution.ParameterFitting(
            PowerLawDistribution.ParameterFitting.Method.MAXIMUM_LIKELIHOOD,
            PowerLawDistribution.DistributionParameters(),
        )
    )


def create_fitting_parameters_stable() -> StableDistribution.FittingParameters:
    """Create fitting parameters for normal distribution."""
    return StableDistribution.FittingParameters(
        StableDistribution.DomainCalculation(),
        StableDistribution.ParameterFitting(
            StableDistribution.ParameterFitting.Method.MLE_SCIPY,
            StableDistribution.DistributionParameters(),
        )
    )


def create_fitting_parameters_uniform() -> UniformDistribution.FittingParameters:
    """Create fitting parameters for normal distribution."""
    return UniformDistribution.FittingParameters(
        UniformDistribution.DomainCalculation(),
        UniformDistribution.ParameterFitting(
            UniformDistribution.ParameterFitting.Method.DEFAULT,
            UniformDistribution.DistributionParameters(),
        )
    )
