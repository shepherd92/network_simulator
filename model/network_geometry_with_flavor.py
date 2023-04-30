#!/usr/bin/env python3
"""This module represents the network geometry with flavor network model.

see: Bianconi: Network geometry with flavor: From complexity to quantum geometry
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from data_set.data_set import DataSet
from model.model import Model
from network.finite_network import FiniteNetwork
from network.property import BaseNetworkProperty


class NetworkGeometryWithFlavorModel(Model):
    """Base class representing a data set."""

    @dataclass
    class Parameters(Model.Parameters):
        """Contain all necessary parameters to construct an ErdosRenyiModel."""

        simplex_dimension: int = 0
        beta: float = 0.1
        flavor: int = -1

    def __init__(self) -> None:
        """Create a network model with default parameters."""
        self._parameters = NetworkGeometryWithFlavorModel.Parameters()

    def set_relevant_parameters_from_data_set(self, data_set: DataSet) -> None:
        """Set the model parameters based ona a data set."""
        num_of_nodes: int = data_set.calc_base_property(BaseNetworkProperty(
            BaseNetworkProperty.Type.NUM_OF_NODES
        ))

        # pylint: disable=attribute-defined-outside-init
        self._parameters.max_dimension = data_set.max_dimension
        self._parameters.num_nodes = num_of_nodes
        # pylint: enable=attribute-defined-outside-init
        self._parameters.simplex_dimension = data_set.max_dimension

    def generate_finite_network(self, _: int | None = None) -> FiniteNetwork:
        """Build a network of the model."""
        assert isinstance(self.parameters, NetworkGeometryWithFlavorModel.Parameters), \
            f'Wrong model parameter type {type(self.parameters)}'
        assert self.parameters.max_dimension >= self.parameters.simplex_dimension, \
            f'The dimension of the simplices ({self.parameters.simplex_dimension}) ' + \
            f'are greater than the dimension that can be handled ({self.parameters.max_dimension}).'

        dimension = self.parameters.simplex_dimension  # alias

        # initialize nodes
        node_ids = np.array(range(self.parameters.num_nodes))
        node_energies = self._generate_node_energies()

        # initialize simplices, all simplices have dimension simplex_dimension
        num_of_simplices = self.parameters.num_nodes - self.parameters.simplex_dimension
        simplices = np.zeros((num_of_simplices, self.parameters.simplex_dimension + 1), dtype=int)
        simplices[0, :] = np.array(node_ids[:self.parameters.simplex_dimension + 1])

        # initialize faces
        num_of_faces = 1 + \
            self.parameters.simplex_dimension * \
            (self.parameters.num_nodes - self.parameters.simplex_dimension)
        # all faces have dimension simplex_dimension - 1
        faces = np.zeros((num_of_faces, self.parameters.simplex_dimension), dtype=int)
        new_faces = np.array([
            face for face in combinations(simplices[0, :], dimension)
        ])
        faces[:dimension + 1, :] = new_faces
        face_energies = np.zeros((num_of_faces))
        energies_of_new_faces = self._calc_simplex_energies(new_faces, node_energies)
        face_energies[:dimension + 1] = energies_of_new_faces

        # initialize n_alpha_values
        # the n_alpha parameter of a face is the number of higher dimensional simplices
        # incident to it - 1
        n_alpha_values = np.zeros(len(faces), dtype=int)

        # initialize face energies
        unnormalized_probabilities_of_new_faces = self._calc_unnormalized_probabilities(
            energies_of_new_faces,
            np.zeros(len(new_faces))
        )
        unnormalized_probabilities = np.zeros(len(faces), dtype=int)
        unnormalized_probabilities[:len(new_faces)] = unnormalized_probabilities_of_new_faces

        # initialize some arrays for easier indexing
        next_simplex_indices = node_ids - dimension
        first_indices_of_new_faces = 1 + dimension * (node_ids - dimension)

        for next_node_id in node_ids[dimension + 1:]:

            probabilities = unnormalized_probabilities / np.sum(unnormalized_probabilities)
            new_faces_first_index = first_indices_of_new_faces[next_node_id]

            # first choose a face with the given probabilities
            chosen_face_index = np.random.choice(range(len(faces)), replace=False, p=probabilities)
            chosen_face = faces[chosen_face_index]

            # add new simplex with the node and the chosen face
            # keep nodes ordered: the index of the new node is always the last
            new_simplex = np.append(chosen_face, next_node_id)
            simplices[next_simplex_indices[next_node_id], :] = new_simplex

            # calculate which new faces will come in
            new_faces = np.array([
                face
                for face in combinations(new_simplex, dimension)
                if next_node_id in face  # the chosen face is not new
            ])
            faces[new_faces_first_index:new_faces_first_index + len(new_faces), :] = new_faces

            # update n_alpha_values (new faces all have n_alpha = 0, do nothing with those)
            n_alpha_values[chosen_face_index] += 1

            # update face_energies
            energies_of_new_faces = self._calc_simplex_energies(new_faces, node_energies)
            face_energies[new_faces_first_index:new_faces_first_index + len(new_faces)] = \
                energies_of_new_faces

            # update unnormalized probabilities
            unnormalized_probabilities[chosen_face_index] = self._calc_unnormalized_probabilities(
                np.expand_dims(face_energies[chosen_face_index], axis=0),
                np.expand_dims(n_alpha_values[chosen_face_index], axis=0)
            )
            unnormalized_probabilities_of_new_faces = self._calc_unnormalized_probabilities(
                energies_of_new_faces,
                np.zeros(len(new_faces))
            )
            unnormalized_probabilities[
                new_faces_first_index:new_faces_first_index + len(new_faces)
            ] = unnormalized_probabilities_of_new_faces

        network = FiniteNetwork(self.parameters.max_dimension)
        network._interactions = simplices
        network._facets = []
        network.add_simplices(simplices)
        network.generate_graph_from_simplicial_complex()

        return network

    def _generate_node_energies(self) -> np.ndarray:
        node_energies = np.random.random_integers(low=0, high=10, size=self.parameters.num_nodes)
        return node_energies

    def _calc_unnormalized_probabilities(
        self,
        simplex_energies: np.ndarray,
        n_alpha_values: np.ndarray
    ) -> np.ndarray:

        unnormalized_probabilities = \
            np.exp(self.parameters.beta * simplex_energies) * \
            (1. + self.parameters.flavor * n_alpha_values)

        return unnormalized_probabilities

    @staticmethod
    def _calc_simplex_energies(simplices: np.ndarray, node_energies: np.ndarray) -> np.ndarray:

        simplex_energies = np.sum(node_energies[simplices], axis=1)
        return simplex_energies

    @property
    def parameters(self) -> NetworkGeometryWithFlavorModel.Parameters:
        """Return the parameters of the network model."""
        return self._parameters

    @parameters.setter
    def parameters(self, value: NetworkGeometryWithFlavorModel.Parameters) -> None:
        self._parameters = value
