#!/usr/bin/env python3
"""Configuration of the Price model."""

from data_set.data_set import DataSet
from model.graph_models.price import PriceModel
from network.property import BaseNetworkProperty
from optimizer.parameter_option import ModelParameterOptions, ParameterOption


PRICE_MODEL_PARAMETERS = PriceModel.Parameters(
    max_dimension=2,
    num_nodes=1000,
    probability_degree_constant=1.0,
)


class PriceParameterOptions(ModelParameterOptions):
    """Represent the parameter options for fitting for Price model."""

    def __init__(self, data_set: DataSet) -> None:
        """Construct parameter options for model fitting."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty.Type.NUM_OF_NODES)

        self.options = (
            #                                        name,           initial guess,  lb,   ub, integer, optimize,  step)
            ParameterOption('max_dimension',  data_set.max_dimension,  0,   10,   True,    False,    1.),
            ParameterOption('num_nodes',            num_of_nodes,  0,    0,   True,    False,    1.),
            ParameterOption('probability_degree_constant',                      1.,  1,  100,   False,    False,    1.),
        )
