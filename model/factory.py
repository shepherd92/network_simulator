#!/usr/bin/env python3
"""This module provides factory methods for creating network models."""

from model.model import Model
from model.graph_models.erdos_renyi import ErdosRenyiModel
from model.graph_models.age_dependent_random_connection import AgeDependentRandomConnectionModel
from model.graph_models.preferential_attachment import PreferentialAttachmentModel
from model.graph_models.price import PriceModel
from model.graph_models.watts_strogatz import WattsStrogatzModel
from model.simplicial_complex_models.age_dependent_random_simplex import AgeDependentRandomSimplexModel
from model.simplicial_complex_models.network_geometry_with_flavor import NetworkGeometryWithFlavorModel

from config_files.model_config.erdos_renyi_config import ERDOS_RENYI_MODEL_PARAMETERS
from config_files.model_config.watts_strogatz_config import WATTS_STROGATZ_PARAMETERS
from config_files.model_config.age_dependent_random_connection_config import (
    AGE_DEPENDENT_RANDOM_CONNECTION_MODEL_PARAMETERS
)
from config_files.model_config.preferential_attachment_config import PREFERENTIAL_ATTACHMENT_MODEL_PARAMETERS
from config_files.model_config.price_config import PRICE_MODEL_PARAMETERS
from config_files.model_config.network_geometry_with_flavor_config import NETWORK_GEOMETRY_WITH_FLAVOR_MODEL_PARAMETERS
from config_files.model_config.age_dependent_random_simplex_config import AGE_DEPENDENT_RANDOM_SIMPLEX_MODEL_PARAMETERS


def create_model(model_type: Model.Type) -> Model:
    """Build the specified model using this factory method."""
    model = Model()
    if model_type == Model.Type.ERDOS_RENYI:
        model = ErdosRenyiModel()
    elif model_type == Model.Type.WATTS_STROGATZ:
        model = WattsStrogatzModel()
    elif model_type == Model.Type.PREFERENTIAL_ATTACHMENT:
        model = PreferentialAttachmentModel()
    elif model_type == Model.Type.PRICE:
        model = PriceModel()
    elif model_type == Model.Type.AGE_DEPENDENT_RANDOM_CONNECTION:
        model = AgeDependentRandomConnectionModel()
    elif model_type == Model.Type.NETWORK_GEOMETRY_WITH_FLAVOR:
        model = NetworkGeometryWithFlavorModel()
    elif model_type == Model.Type.AGE_DEPENDENT_RANDOM_SIMPLEX:
        model = AgeDependentRandomSimplexModel()
    else:
        raise NotImplementedError(f'The requested model type {model_type.name} is not implemented.')

    return model


def load_default_parameters(model_type: Model.Type) -> Model.Parameters:
    """Load parameters for the specified model."""
    parameters: Model.Parameters = Model.Parameters()
    if model_type == Model.Type.ERDOS_RENYI:
        parameters = ERDOS_RENYI_MODEL_PARAMETERS
    elif model_type == Model.Type.WATTS_STROGATZ:
        parameters = WATTS_STROGATZ_PARAMETERS
    elif model_type == Model.Type.PREFERENTIAL_ATTACHMENT:
        parameters = PREFERENTIAL_ATTACHMENT_MODEL_PARAMETERS
    elif model_type == Model.Type.PRICE:
        parameters = PRICE_MODEL_PARAMETERS
    elif model_type == Model.Type.AGE_DEPENDENT_RANDOM_CONNECTION:
        parameters = AGE_DEPENDENT_RANDOM_CONNECTION_MODEL_PARAMETERS
    elif model_type == Model.Type.NETWORK_GEOMETRY_WITH_FLAVOR:
        parameters = NETWORK_GEOMETRY_WITH_FLAVOR_MODEL_PARAMETERS
    elif model_type == Model.Type.AGE_DEPENDENT_RANDOM_SIMPLEX:
        parameters = AGE_DEPENDENT_RANDOM_SIMPLEX_MODEL_PARAMETERS
    else:
        raise NotImplementedError(
            f'The requested model type {model_type} is not implemented.'
        )

    return parameters
