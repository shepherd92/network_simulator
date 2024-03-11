#!/usr/bin/env python3
"""Configuration of the Age Dependent Random Simplex model."""

import numpy as np

from data_set.data_set import DataSet
from distribution.approximation import DistributionApproximation
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.factory import (
    create_fitting_parameters_normal,
    create_fitting_parameters_power_law_model,
)
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution
from distribution.theoretical.power_law_distribution import PowerLawDistribution
from network.property import BaseNetworkProperty, DerivedNetworkProperty
from optimizer.parameter_option import ModelParameterOptions, ParameterOption


SCALAR_PROPERTY_PARAMS_TO_FIT: tuple[DerivedNetworkProperty, ...] = (
    DerivedNetworkProperty(
        name='num_of_edges',
        source_base_property=BaseNetworkProperty(BaseNetworkProperty.Type.NUM_OF_EDGES),
        fitting_parameters=create_fitting_parameters_normal(),
        theoretical_approximation_type=TheoreticalDistribution.Type.NORMAL,
    ),
)


class AgeDependentRandomSimplexParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Age-Dependent Random Simplex model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        super().__init__()
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_VERTICES)
        average_degree: float = data_set.calc_base_property(BaseNetworkProperty.Type.AVERAGE_DEGREE)
        in_degree_distribution: EmpiricalDistribution = data_set.calc_base_property(
            BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION
        )
        approximation = DistributionApproximation(
            in_degree_distribution,
            TheoreticalDistribution.Type.POWER_LAW
        )
        assert isinstance(approximation.theoretical, PowerLawDistribution)
        approximation.fit(create_fitting_parameters_power_law_model())
        gamma_guess = 1. / (approximation.theoretical.parameters.exponent - 1.)
        beta_guess = (1. - gamma_guess) * average_degree

        self.options = (
            #               name,               initial guess,           lb,         ub,  integer, optimize,  step)
            ParameterOption('max_dimension',    data_set.max_dimension,  0,          10,     True,    False,    1.),
            ParameterOption('num_nodes',        num_of_nodes,            0,           0,     True,    False,    1.),
            ParameterOption('torus_dimension',  1,                       1,           5,     True,    False,    1.),
            ParameterOption('alpha',            0.5,                     0.5,  np.infty,    False,    False,    1e-3),
            ParameterOption('beta',             beta_guess,              0.,   np.infty,    False,     True,    1e-1),
            ParameterOption('gamma',            gamma_guess,             0.,        1.0,    False,    False,    1e-3),
        )


class ErdosRenyiParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Erdos-Renyi model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        super().__init__()
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_VERTICES)
        num_of_edges: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_EDGES)
        edge_probability_guess = num_of_edges / (num_of_nodes * (num_of_nodes - 1) / 2)

        self.options = (
            #                name,              initial guess,           lb,       ub,   integer, optimize, step)
            ParameterOption('max_dimension',    data_set.max_dimension,  0,        10,   True,    False,    1.),
            ParameterOption('num_nodes',        num_of_nodes,            0,         0,   True,    False,    1.),
            ParameterOption('edge_probability', edge_probability_guess,  0.0,       1.0, False,   False,    1.),
        )


class NetworkGeometryWithFlavorParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Network Geometry With Flavor model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        super().__init__()
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_VERTICES)

        self.options = (
            #                              name,           initial guess,  lb,    ub, integer, optimize,  step)
            ParameterOption('max_dimension',      data_set.max_dimension,  0,    10,   True,    False,  1.),
            ParameterOption('num_nodes',          num_of_nodes,            0,     0,   True,    False,  1.),
            ParameterOption('simplex_dimension',  data_set.max_dimension,  0,     5,   True,    False,  1.),
            ParameterOption('beta',               1.,                      0.,  100.,  False,   False,  1e-2),
            ParameterOption('flavor',             0,                      -1,     5,   True,    False,  1.),
        )


class PreferentialAttachmentParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Preferential Attachment model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        super().__init__()
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_VERTICES)
        num_of_edges: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_EDGES)
        edges_of_new_node_guess = round(num_of_edges / num_of_nodes)

        self.options = (
            #                              name,            initial guess,  lb,   ub, integer, optimize,  step)
            ParameterOption('max_dimension',   data_set.max_dimension,  0,   10,   True,    False,  1.),
            ParameterOption('num_nodes',             num_of_nodes,  0,    0,   True,    False,  1.),
            ParameterOption('edges_of_new_node',  edges_of_new_node_guess,  1,  100,   True,    False,  1.),
        )


class PriceParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Price model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        super().__init__()
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_VERTICES)

        self.options = (
            #                                        name,           initial guess,  lb,   ub, integer,  optimize, step)
            ParameterOption('max_dimension',               data_set.max_dimension,  0,   10,   True,     False,     1.),
            ParameterOption('num_nodes',                   num_of_nodes,            0,    0,   True,     False,     1.),
            ParameterOption('probability_degree_constant', 1.,                      1,  100,   False,    False,     1.),
        )


class WattsStrogatzParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Watts-Strogatz model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        super().__init__()
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_VERTICES)
        num_of_edges: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_EDGES)
        edges_of_new_node_guess = round(num_of_edges / num_of_nodes)

        self.options = (
            #                                 name,            initial guess,  lb,  ub, integer, optimize,  step)
            ParameterOption('max_dimension',   data_set.max_dimension,  0,  10,   True,    False,  1.),
            ParameterOption('num_nodes',             num_of_nodes,  0,   0,   True,    False,  1.),
            ParameterOption('edges_of_new_node',  edges_of_new_node_guess,  0,   5,   True,    False,  1.),
            ParameterOption('rewiring_probability',                      0.2,  0.,  1.,   False,    False,  1e-2),
        )
