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
        properties_to_calculate_finite: list[BaseNetworkProperty] = [
            BaseNetworkProperty.num_of_vertices,
            BaseNetworkProperty.num_of_edges,
            BaseNetworkProperty.num_of_triangles,
            # BaseNetworkProperty.num_of_interactions,
            # BaseNetworkProperty.edges,
            BaseNetworkProperty.mean_degree,
            BaseNetworkProperty.max_degree,
            # BaseNetworkProperty.mean_clustering,
            BaseNetworkProperty.num_of_connected_components,
            BaseNetworkProperty.interaction_vertex_degree_distribution,
            BaseNetworkProperty.simplex_dimension_distribution,
            # BaseNetworkProperty.facet_dimension_distribution,
            BaseNetworkProperty.vertex_edge_degree_distribution,
            # BaseNetworkProperty.in_degree_distribution,
            # BaseNetworkProperty.out_degree_distribution,
            BaseNetworkProperty.vertex_interaction_degree_distribution,
            BaseNetworkProperty.edge_triangle_degree_distribution,
            # BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
            BaseNetworkProperty.betti_numbers,
            # BaseNetworkProperty.betti_numbers_by_component,
            BaseNetworkProperty.num_of_vertices_by_component,
            # BaseNetworkProperty.persistence,
            # BaseNetworkProperty.persistence_pairs,
        ]
        properties_to_calculate_infinite: list[BaseNetworkProperty] = [
            # BaseNetworkProperty.vertex_edge_degree_distribution,
            # BaseNetworkProperty.in_degree_distribution,
            # BaseNetworkProperty.out_degree_distribution,
            # BaseNetworkProperty.edge_triangle_degree_distribution,
            # BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
        ]

    class Fitting(NamedTuple):
        """Model fitting configuration."""

    class Testing(NamedTuple):
        """Model testing configuration for networks."""

        mode: Model.Mode = Model.Mode.INFINITE
        num_of_simulations: int = 10
        num_of_infinite_networks: int = 10

    type_: Model.Type = Model.Type.AGE_DEPENDENT_RANDOM_HYPERGRAPH
    set_params_from_data_set: bool = False
    analysis = Analysis()
    fitting = Fitting()
    network_testing = Testing()


# ==============================================================================
NETWORK_MAGNITUDE = 4
GAMMA = 0.7
GAMMA_PRIME = 0.2
EXPECTED_VERTEX_INTERACTION_DEGREE = 3
HYPERGRAPH_MODEL_PARAMETERS = HypergraphModel.Parameters(
    max_dimension=2,
    network_size=10**NETWORK_MAGNITUDE,  # vertex_intensity: expected number of vertices
    interaction_intensity=10**NETWORK_MAGNITUDE,  # interaction_intensity: expected number of interactions
    beta=EXPECTED_VERTEX_INTERACTION_DEGREE * 0.5 * (1. - GAMMA) * (1. - GAMMA_PRIME) * 10**(-NETWORK_MAGNITUDE),
    # beta=10**(-NETWORK_MAGNITUDE),
    # network_size=,
    # interaction_intensity=,
    # beta=,
    gamma=GAMMA,
    gamma_prime=GAMMA_PRIME,
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
