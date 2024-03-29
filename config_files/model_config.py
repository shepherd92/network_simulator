#!/usr/bin/env python3
"""Configuration of the Age Dependent Random Simplex model."""

from dataclasses import dataclass
from typing import NamedTuple

from model.model import Model
from model.age_dependent_random_simplex import AgeDependentRandomSimplexModel
from model.hypergraph import HypergraphModel
from model.erdos_renyi import ErdosRenyiModel
from model.network_geometry_with_flavor import NetworkGeometryWithFlavorModel
from model.preferential_attachment import PreferentialAttachmentModel
from model.price import PriceModel
from model.watts_strogatz import WattsStrogatzModel
from network.property import BaseNetworkProperty


@dataclass
class ModelConfig:
    """Model configuration."""

    class Analysis(NamedTuple):
        """Model analysis configuration."""

        component_index_from_largest: int = -1
        plot: bool = True
        num_of_infinite_networks: int = 0
        properties_to_calculate_finite: list[BaseNetworkProperty.Type] = [
            BaseNetworkProperty.Type.NUM_OF_VERTICES,
            BaseNetworkProperty.Type.NUM_OF_EDGES,
            BaseNetworkProperty.Type.NUM_OF_TRIANGLES,
            BaseNetworkProperty.Type.NUM_OF_INTERACTIONS,
            # BaseNetworkProperty.Type.EDGES,
            BaseNetworkProperty.Type.AVERAGE_DEGREE,
            BaseNetworkProperty.Type.MAX_DEGREE,
            BaseNetworkProperty.Type.AVG_CLUSTERING,
            BaseNetworkProperty.Type.NUM_OF_CONNECTED_COMPONENTS,
            BaseNetworkProperty.Type.INTERACTION_DIMENSION_DISTRIBUTION,
            # BaseNetworkProperty.Type.SIMPLEX_DIMENSION_DISTRIBUTION,
            # BaseNetworkProperty.Type.FACET_DIMENSION_DISTRIBUTION,
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION,
            BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION,
            BaseNetworkProperty.Type.OUT_DEGREE_DISTRIBUTION,
            BaseNetworkProperty.Type.VERTEX_INTERACTION_DEGREE_DISTRIBUTION,
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1,
            BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2,
            BaseNetworkProperty.Type.BETTI_NUMBERS,
            BaseNetworkProperty.Type.BETTI_NUMBERS_BY_COMPONENT,
            BaseNetworkProperty.Type.VERTICES_BY_COMPONENT,
            # BaseNetworkProperty.Type.PERSISTENCE_PAIRS,
        ]
        properties_to_calculate_infinite: list[BaseNetworkProperty.Type] = [
            BaseNetworkProperty.Type.DEGREE_DISTRIBUTION,
            # BaseNetworkProperty.Type.IN_DEGREE_DISTRIBUTION,
            # BaseNetworkProperty.Type.OUT_DEGREE_DISTRIBUTION,
            # BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_1,
            # BaseNetworkProperty.Type.HIGHER_ORDER_DEGREE_DISTRIBUTION_2,
        ]

    class Fitting(NamedTuple):
        """Model fitting configuration."""

    class Testing(NamedTuple):
        """Model testing configuration for networks."""

        num_of_simulations: int = 100
        num_of_infinite_networks: int = 0

    type_: Model.Type = Model.Type.AGE_DEPENDENT_RANDOM_HYPERGRAPH
    set_params_from_data_set: bool = False
    analysis = Analysis()
    fitting = Fitting()
    network_testing = Testing()


# ==============================================================================
HYPERGRAPH_MODEL_PARAMETERS = HypergraphModel.Parameters(
    max_dimension=2,
    network_size=1e6,  # vertex_intensity: expected number of vertices
    interaction_intensity=1e6,  # interaction_intensity: expected number of interactions
    beta=1e-6,
    gamma=0.7,
    gamma_prime=0.2,
)
# ==============================================================================


AGE_DEPENDENT_RANDOM_SIMPLEX_MODEL_PARAMETERS = AgeDependentRandomSimplexModel.Parameters(
    max_dimension=2,
    network_size=100,  # expected number of nodes
    alpha=0.5,
    beta=1.0,
    gamma=0.7,
)


ERDOS_RENYI_MODEL_PARAMETERS = ErdosRenyiModel.Parameters(
    max_dimension=2,
    network_size=1000,
    edge_probability=0.5,
)


NETWORK_GEOMETRY_WITH_FLAVOR_MODEL_PARAMETERS = NetworkGeometryWithFlavorModel.Parameters(
    max_dimension=2,
    network_size=1000,
    simplex_dimension=2,
    beta=1.,
    flavor=0,
)


PREFERENTIAL_ATTACHMENT_MODEL_PARAMETERS = PreferentialAttachmentModel.Parameters(
    max_dimension=2,
    network_size=1000,
    edges_of_new_node=2,
)


PRICE_MODEL_PARAMETERS = PriceModel.Parameters(
    max_dimension=2,
    network_size=1000,
    probability_degree_constant=1.0,
)


WATTS_STROGATZ_PARAMETERS = WattsStrogatzModel.Parameters(
    max_dimension=2,
    network_size=1000,
    edges_of_new_node=2,
    rewiring_probability=0.2,
)
