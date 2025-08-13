#!/usr/bin/env python3
"""This module provides factory methods for creating network models."""

from model.model import Model
from model.erdos_renyi import ErdosRenyiModel
from model.adrcm import AdrcmModel
from model.hypergraph import HypergraphModel

from config_files.model_config import ERDOS_RENYI_MODEL_PARAMETERS
from config_files.model_config import ADRCM_MODEL_PARAMETERS
from config_files.model_config import HYPERGRAPH_MODEL_PARAMETERS


def create_model(model_type: Model.Type) -> Model:
    """Build the specified model using this factory method."""
    model = Model()
    if model_type == Model.Type.ERDOS_RENYI:
        model = ErdosRenyiModel()
    elif model_type == Model.Type.ADRCM:
        model = AdrcmModel()
    elif model_type == Model.Type.HYPERGRAPH:
        model = HypergraphModel()
    else:
        raise NotImplementedError(f'The requested model type {model_type.name} is not implemented.')

    return model


def load_default_parameters(model_type: Model.Type) -> Model.Parameters:
    """Load parameters for the specified model."""
    parameters: Model.Parameters = Model.Parameters()
    if model_type == Model.Type.ERDOS_RENYI:
        parameters = ERDOS_RENYI_MODEL_PARAMETERS
    elif model_type == Model.Type.ADRCM:
        parameters = ADRCM_MODEL_PARAMETERS
    elif model_type == Model.Type.HYPERGRAPH:
        parameters = HYPERGRAPH_MODEL_PARAMETERS
    else:
        raise NotImplementedError(
            f'The requested model type {model_type} is not implemented.'
        )

    return parameters
