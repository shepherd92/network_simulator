#!/usr/bin/env python3
"""Configuration of the Age Dependent Random Simplex model."""

from dataclasses import dataclass
from typing import NamedTuple

from model.model import Model
from model.adrcm import AdrcmModel
from model.hypergraph import HypergraphModel
from model.erdos_renyi import ErdosRenyiModel
from network.property import BaseNetworkProperty


@dataclass
class ModelConfig:
    """Model configuration."""

    class Analysis(NamedTuple):
        """Model analysis configuration."""

        class Finite(NamedTuple):
            """Finite network analysis configuration."""
            enable: bool = True
            component_index_from_largest: int = -1

            properties_to_calculate: list[BaseNetworkProperty] = [
                BaseNetworkProperty.num_of_vertices,
                # BaseNetworkProperty.num_of_edges,
                BaseNetworkProperty.num_of_triangles,
                BaseNetworkProperty.num_of_interactions,
                BaseNetworkProperty.mean_degree,
                BaseNetworkProperty.max_degree,
                BaseNetworkProperty.interaction_vertex_degree_distribution,
                BaseNetworkProperty.simplex_dimension_distribution,
                BaseNetworkProperty.vertex_edge_degree_distribution,
                BaseNetworkProperty.vertex_interaction_degree_distribution,
                BaseNetworkProperty.edge_interaction_degree_distribution,
                BaseNetworkProperty.triangle_interaction_degree_distribution,
                BaseNetworkProperty.edge_triangle_degree_distribution,
                BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
                BaseNetworkProperty.betti_numbers,
                BaseNetworkProperty.betti_numbers_by_component,
                BaseNetworkProperty.num_of_vertices_by_component,
                # BaseNetworkProperty.persistence_intervals,
                # BaseNetworkProperty.persistence_pairs,
            ]

        class InfiniteSet(NamedTuple):
            """Infinite network set analysis configuration."""
            enable: bool = False
            num_of_infinite_networks: int = 100
            properties_to_calculate: list[BaseNetworkProperty] = [
                # BaseNetworkProperty.vertex_edge_degree_distribution,
                # BaseNetworkProperty.edge_triangle_degree_distribution,
                # BaseNetworkProperty.vertex_interaction_degree_distribution,
                # BaseNetworkProperty.edge_interaction_degree_distribution,
                # BaseNetworkProperty.triangle_interaction_degree_distribution,
                # BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
            ]

        class Infinite(NamedTuple):
            """Infinite network analysis configuration."""
            enable: bool = True
            plot: bool = True
            typical_mark: float = 0.1

        finite: Finite = Finite()
        infinite: Infinite = Infinite()
        infinite_set: InfiniteSet = InfiniteSet()

        power_law_fitting_minimum_value: float = 10

        # plotting
        plot: bool = True
        plot_entire_network: bool = True
        plot_network_giant_component: bool = True
        plot_network_determined_positions: bool = True

    class Testing(NamedTuple):
        """Model testing configuration for networks."""

        mode: Model.Mode = Model.Mode.FINITE
        num_of_simulations: int = 100
        num_of_infinite_networks: int = 1000

    type_: Model.Type = Model.Type.HYPERGRAPH
    set_params_from_data_set: bool = False
    analysis = Analysis()
    network_testing = Testing()


# ==============================================================================
NETWORK_MAGNITUDE = 2
GAMMA = 0.7
GAMMA_PRIME = 0.2
EXPECTED_VERTEX_INTERACTION_DEGREE = 3.

HYPERGRAPH_MODEL_PARAMETERS = HypergraphModel.Parameters(
    max_dimension=2,
    network_size=10**NETWORK_MAGNITUDE,  # vertex_intensity: expected number of vertices
    interaction_intensity=10**NETWORK_MAGNITUDE,  # interaction_intensity: expected number of interactions
    beta=EXPECTED_VERTEX_INTERACTION_DEGREE * 0.5 * (1. - GAMMA) * (1. - GAMMA_PRIME) * 10**(-NETWORK_MAGNITUDE),
    gamma=GAMMA,
    gamma_prime=GAMMA_PRIME,
    weighted=False,
    interactions_enough=False,
)
# ==============================================================================


ADRCM_MODEL_PARAMETERS = AdrcmModel.Parameters(
    max_dimension=2,
    network_size=10**NETWORK_MAGNITUDE,  # expected number of nodes
    alpha=0.5,
    beta=1.0,
    gamma=0.7,
)


ERDOS_RENYI_MODEL_PARAMETERS = ErdosRenyiModel.Parameters(
    max_dimension=2,
    network_size=1000,
    edge_probability=0.5,
)
