#!/usr/bin/env python3
"""Unittest for InfiniteHypergraph."""

import unittest

from network.infinite_hypergraph import InfiniteHypergraph
from network.infinite_network import CppInfiniteHypergraph
from network.property import BaseNetworkProperty

class InfiniteHypergraphTest(unittest.TestCase):
    """Test for the finite clique complex."""

    def setUp(self) -> None:
        cpp_network = CppInfiniteHypergraph(
            3, # max_dimension
            [0, 1, 2, 3, 4, 5, 6], # vertices
            [
                [0, 1, 2, 3],
                [0, 2, 3],
                [4, 5],
            ], # interactions
            0.45, # typical_mark
            [0.7, 0.5, 0.5, 0.6, 0.8, 0.3, 0.8], # mark_list
        )
        self.network = InfiniteHypergraph(cpp_network)

    def test_number_of_simplices(self) -> None:
        # typical vertex is not included in the simplices
        self.assertEqual(self.network.num_simplices(0), 0)
        self.assertEqual(self.network.num_simplices(1), 6)
        self.assertEqual(self.network.num_simplices(2), 6)
        self.assertEqual(self.network.num_simplices(3), 4)

    def test_degree_distribution_0_1(self) -> None:
        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.vertex_edge_degree_distribution)
        self.assertEqual(sorted(degrees), [7])

    def test_degree_distribution_1_2(self) -> None:
        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.edge_triangle_degree_distribution)
        self.assertEqual(sorted(degrees), [0, 1, 3, 3, 3, 3])

    def test_degree_distribution_2_3(self) -> None:
        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.triangle_tetrahedra_degree_distribution)
        self.assertEqual(sorted(degrees), [2, 2, 2, 2, 2, 2])

    def test_typical_vertex_interaction_degree(self) -> None:
        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.vertex_interaction_degree_distribution)
        self.assertEqual(sorted(degrees), [3])

    def test_typical_edge_interaction_degree(self) -> None:
        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.edge_interaction_degree_distribution)
        self.assertEqual(sorted(degrees), [0, 1, 1, 2, 2, 2])

    def test_typical_triangle_interaction_degree(self) -> None:
        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.triangle_interaction_degree_distribution)
        self.assertEqual(sorted(degrees), [1, 1, 1, 2, 2, 2])

    def test_interaction_vertex_degree_distribution(self) -> None:
        distribution: list[int] = self.network.calc_base_property(BaseNetworkProperty.interaction_vertex_degree_distribution)
        self.assertEqual(sorted(distribution), [2, 3, 4])


if __name__ == '__main__':
    unittest.main()
