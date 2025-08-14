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
                BaseNetworkProperty.num_of_edges,
                BaseNetworkProperty.num_of_triangles,
                BaseNetworkProperty.mean_degree,
                BaseNetworkProperty.max_degree,
                BaseNetworkProperty.simplex_dimension_distribution,
                BaseNetworkProperty.vertex_edge_degree_distribution,
                BaseNetworkProperty.edge_triangle_degree_distribution,
                BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
                BaseNetworkProperty.betti_numbers,
                BaseNetworkProperty.betti_numbers_by_component,
                BaseNetworkProperty.num_of_vertices_by_component,
                # Below properties are only available for hypergraphs
                BaseNetworkProperty.num_of_interactions,
                BaseNetworkProperty.interaction_vertex_degree_distribution,
                BaseNetworkProperty.vertex_interaction_degree_distribution,
                BaseNetworkProperty.edge_interaction_degree_distribution,
                BaseNetworkProperty.triangle_interaction_degree_distribution,
                BaseNetworkProperty.persistence_intervals,
                BaseNetworkProperty.persistence_pairs,
            ]

        class Infinite(NamedTuple):
            """Infinite network analysis configuration."""
            enable: bool = False
            num_of_infinite_networks: int = 100
            properties_to_calculate: list[BaseNetworkProperty] = [
                BaseNetworkProperty.vertex_edge_degree_distribution,
                BaseNetworkProperty.edge_triangle_degree_distribution,
                BaseNetworkProperty.triangle_tetrahedra_degree_distribution,
                # Below properties are only available for hypergraphs
                BaseNetworkProperty.vertex_interaction_degree_distribution,
                BaseNetworkProperty.edge_interaction_degree_distribution,
                BaseNetworkProperty.triangle_interaction_degree_distribution,
            ]

        finite: Finite = Finite()
        infinite: Infinite = Infinite()

        # plotting
        plot: bool = True
        plot_entire_network: bool = False
        plot_network_giant_component: bool = False
        plot_network_determined_positions: bool = True

    class Testing(NamedTuple):
        """Model testing configuration for networks."""

        mode: Model.Mode = Model.Mode.FINITE
        num_of_simulations: int = 100
        num_of_infinite_networks: int = 100

    type_: Model.Type = Model.Type.HYPERGRAPH
    set_params_from_data_set: bool = False
    analysis = Analysis()
    network_testing = Testing()

    # ==============================================================================
    max_dimension: int = 2
    network_magnitude = 6
    gamma = 0.7
    gamma_prime = 0.2
    expected_vertex_interaction_degree = 3.
    # ==============================================================================
    # ==============================================================================
    hypergraph_model_parameters = HypergraphModel.Parameters(
        max_dimension=max_dimension,
        network_size=10**network_magnitude,  # vertex_intensity: expected number of vertices
        interaction_intensity=10**network_magnitude,  # interaction_intensity: expected number of interactions
        beta=expected_vertex_interaction_degree * 0.5 * (1. - gamma) * (1. - gamma_prime) * 10**(-network_magnitude),
        gamma=gamma,
        gamma_prime=gamma_prime,
        weighted=False,
        interactions_enough=False,
    )
    # ==============================================================================
    adrcm_model_parameters = AdrcmModel.Parameters(
        max_dimension=max_dimension,
        network_size=10**network_magnitude,  # expected number of nodes
        alpha=0.5,
        beta=1.0,
        gamma=gamma,
    )
    # ==============================================================================
    erdos_renyi_model_parameters = ErdosRenyiModel.Parameters(
        max_dimension=max_dimension,
        network_size=10**network_magnitude,
        edge_probability=0.5,
    )
