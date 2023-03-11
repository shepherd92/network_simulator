#!/usr/bin/env python3

from model.model import Model

from config_files.model_config.erdos_renyi_config import ErdosRenyiParameterOptions
from config_files.model_config.watts_strogatz_config import WattsStrogatzParameterOptions
from config_files.model_config.age_dependent_random_connection_config import (
    AgeDependentRandomConnectionParameterOptions
)
from config_files.model_config.preferential_attachment_config import PreferentialAttachmentParameterOptions
from config_files.model_config.price_config import PriceParameterOptions
from config_files.model_config.network_geometry_with_flavor_config import NetworkGeometryWithFlavorParameterOptions
from config_files.model_config.age_dependent_random_simplex_config import AgeDependentRandomSimplexParameterOptions
from data_set.data_set import DataSet
from optimizer.parameter_option import ModelParameterOptions


def create_parameter_options(model_type: Model.Type, data_set: DataSet) -> ModelParameterOptions:
    """Build parameter options for the specified model."""

    parameter_options: ModelParameterOptions = ModelParameterOptions()
    if model_type == Model.Type.ERDOS_RENYI:
        parameter_options = ErdosRenyiParameterOptions(data_set)
    elif model_type == Model.Type.WATTS_STROGATZ:
        parameter_options = WattsStrogatzParameterOptions(data_set)
    elif model_type == Model.Type.PREFERENTIAL_ATTACHMENT:
        parameter_options = PreferentialAttachmentParameterOptions(data_set)
    elif model_type == Model.Type.PRICE:
        parameter_options = PriceParameterOptions(data_set)
    elif model_type == Model.Type.AGE_DEPENDENT_RANDOM_CONNECTION:
        parameter_options = AgeDependentRandomConnectionParameterOptions(data_set)
    elif model_type == Model.Type.NETWORK_GEOMETRY_WITH_FLAVOR:
        parameter_options = NetworkGeometryWithFlavorParameterOptions(data_set)
    elif model_type == Model.Type.AGE_DEPENDENT_RANDOM_SIMPLEX:
        parameter_options = AgeDependentRandomSimplexParameterOptions(data_set)
    else:
        raise NotImplementedError(
            f'The requested model type {model_type} is not implemented.'
        )

    return parameter_options
