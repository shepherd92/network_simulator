#!/usr/bin/env python3
"""Configuration of the Age Dependent Random Connection model."""

import numpy as np

from data_set.data_set import DataSet
from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from model.graph_models.age_dependent_random_connection import AgeDependentRandomConnectionModel
from network.property import BaseNetworkProperty
from optimizer.parameter_option import ModelParameterOptions, ParameterOption


AGE_DEPENDENT_RANDOM_CONNECTION_MODEL_PARAMETERS = AgeDependentRandomConnectionModel.Parameters(
    max_dimension=2,
    num_nodes=1000,
    torus_dimension=1,
    alpha=0.5,
    beta=1.0,
    gamma=0.5,
)


class AgeDependentRandomConnectionParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Age-Dependent Random Connection model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)
        edge_density: float = data_set.calc_base_property(BaseNetworkProperty.Type.AVERAGE_DEGREE)
        in_degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION
        )
        approximation = DistributionApproximation(in_degree_distribution, TheoreticalDistribution.Type.POWER_LAW)
        assert isinstance(approximation.theoretical, PowerLawDistribution)
        approximation.fit()
        gamma_guess = 1. / (approximation.theoretical.parameters.exponent - 1.)
        beta_guess = (1. - gamma_guess) * edge_density

        self.options = (
            #               name,              initial guess,          lb,          ub,  integer,  optimize,  step)
            ParameterOption('max_dimension',   data_set.max_dimension,  0,          10,     True,     False,    1.),
            ParameterOption('num_nodes',       num_of_nodes,            0,           0,     True,     False,    1.),
            ParameterOption('torus_dimension', 1,                       1,           5,     True,     False,    1.),
            ParameterOption('alpha',           0.5,                     0.5,  np.infty,    False,     False,  1e-3),
            ParameterOption('beta',            beta_guess,              0.,   np.infty,    False,      True,  1e-1),
            ParameterOption('gamma',           gamma_guess,             0.,        1.0,    False,     False,  1e-3),
        )
