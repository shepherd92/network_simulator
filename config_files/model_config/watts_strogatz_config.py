#!/usr/bin/env python3
"""Configuration of the Watts Strogatz model."""

from data_set.data_set import DataSet
from model.watts_strogatz import WattsStrogatzModel
from network.property import BaseNetworkProperty
from optimizer.parameter_option import ModelParameterOptions, ParameterOption


WATTS_STROGATZ_PARAMETERS = WattsStrogatzModel.Parameters(
    max_dimension=2,
    num_nodes=1000,
    edges_of_new_node=2,
    rewiring_probability=0.2,
)


class WattsStrogatzParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Watts-Strogatz model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)
        num_of_edges: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_EDGES)
        edges_of_new_node_guess = round(num_of_edges / num_of_nodes)

        self.options = (
            #                                 name,            initial guess,  lb,  ub, integer, optimize,  step)
            ParameterOption('max_dimension',   data_set.max_dimension,  0,  10,   True,    False,  1.),
            ParameterOption('num_nodes',             num_of_nodes,  0,   0,   True,    False,  1.),
            ParameterOption('edges_of_new_node',  edges_of_new_node_guess,  0,   5,   True,    False,  1.),
            ParameterOption('rewiring_probability',                      0.2,  0.,  1.,   False,    False,  1e-2),
        )
