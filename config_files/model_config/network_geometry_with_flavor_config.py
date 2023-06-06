#!/usr/bin/env python3
"""Configuration of the Network Geometry With Flavor model."""

from data_set.data_set import DataSet
from model.network_geometry_with_flavor import NetworkGeometryWithFlavorModel
from network.property import BaseNetworkProperty
from optimizer.parameter_option import ModelParameterOptions, ParameterOption


NETWORK_GEOMETRY_WITH_FLAVOR_MODEL_PARAMETERS = NetworkGeometryWithFlavorModel.Parameters(
    max_dimension=2,
    num_nodes=1000,
    simplex_dimension=2,
    beta=1.,
    flavor=0,
)


class NetworkGeometryWithFlavorParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Network Geometry With Flavor model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)

        self.options = (
            #                              name,           initial guess,  lb,    ub, integer, optimize,  step)
            ParameterOption('max_dimension',      data_set.max_dimension,  0,    10,   True,    False,  1.),
            ParameterOption('num_nodes',          num_of_nodes,            0,     0,   True,    False,  1.),
            ParameterOption('simplex_dimension',  data_set.max_dimension,  0,     5,   True,    False,  1.),
            ParameterOption('beta',               1.,                      0.,  100.,  False,   False,  1e-2),
            ParameterOption('flavor',             0,                      -1,     5,   True,    False,  1.),
        )
