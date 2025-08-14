#!/usr/bin/env python3
"""This module provides factory methods for creating network models."""

from model.model import Model
from model.erdos_renyi import ErdosRenyiModel
from model.adrcm import AdrcmModel
from model.hypergraph import HypergraphModel

from config_files.model_config import ModelConfig


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
        parameters = ModelConfig.erdos_renyi_model_parameters
    elif model_type == Model.Type.ADRCM:
        parameters = ModelConfig.adrcm_model_parameters
    elif model_type == Model.Type.HYPERGRAPH:
        parameters = ModelConfig.hypergraph_model_parameters
    else:
        raise NotImplementedError(
            f'The requested model type {model_type} is not implemented.'
        )

    return parameters
