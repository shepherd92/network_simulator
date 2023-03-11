#!/usr/bin/env python3
"""Configuration of the Erdos - Renyi model."""

from data_set.data_set import DataSet
from model.graph_models.erdos_renyi import ErdosRenyiModel
from network.property import BaseNetworkProperty
from optimizer.parameter_option import ModelParameterOptions, ParameterOption


ERDOS_RENYI_MODEL_PARAMETERS = ErdosRenyiModel.Parameters(
    max_dimension=2,
    num_nodes=1000,
    edge_probability=0.5,
)


class ErdosRenyiParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Erdos-Renyi model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)
        num_of_edges: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_EDGES)
        edge_probability_guess = num_of_edges / (num_of_nodes * (num_of_nodes - 1) / 2)

        self.options = (
            #                             name,           initial guess,   lb,        ub, integer, optimize,  step)
            ParameterOption('max_dimension',  data_set.max_dimension,  0,        10,   True,    False,  1.),
            ParameterOption('num_nodes',            num_of_nodes,  0,         0,   True,    False,  1.),
            ParameterOption('edge_probability',  edge_probability_guess,  0.0,       1.0,   False,    False,  1.),
        )
