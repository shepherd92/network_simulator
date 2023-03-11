#!/usr/bin/env python3
"""Factory methods for creating theoretical distributions."""

from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.normal_distribution import NormalDistribution
from distribution.theoretical.poisson_distribution import PoissonDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from distribution.theoretical.stable_distribution import StableDistribution


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
    raise NotImplementedError
