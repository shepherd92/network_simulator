#!/usr/bin/env python3
"""Factory methods for the model optimizer."""

from model.model import Model

from config_files.model_fitting import AgeDependentRandomSimplexParameterOptions
from config_files.model_fitting import ErdosRenyiParameterOptions
from config_files.model_fitting import NetworkGeometryWithFlavorParameterOptions
from config_files.model_fitting import PreferentialAttachmentParameterOptions
from data_set.data_set import DataSet
from optimizer.parameter_option import ModelParameterOptions


def create_parameter_options(model_type: Model.Type, data_set: DataSet) -> ModelParameterOptions:
    """Build parameter options for the specified model."""
    parameter_options: ModelParameterOptions = ModelParameterOptions()
    if model_type == Model.Type.ERDOS_RENYI:
        parameter_options = ErdosRenyiParameterOptions(data_set)
    elif model_type == Model.Type.PREFERENTIAL_ATTACHMENT:
        parameter_options = PreferentialAttachmentParameterOptions(data_set)
    elif model_type == Model.Type.NETWORK_GEOMETRY_WITH_FLAVOR:
        parameter_options = NetworkGeometryWithFlavorParameterOptions(data_set)
    elif model_type == Model.Type.AGE_DEPENDENT_RANDOM_SIMPLEX:
        parameter_options = AgeDependentRandomSimplexParameterOptions(data_set)
    else:
        raise NotImplementedError(
            f'The requested model type {model_type} is not implemented.'
        )

    return parameter_options
