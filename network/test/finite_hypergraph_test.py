#!/usr/bin/env python3
"""Unittest for FiniteHypergraph."""

import unittest

from distribution.empirical_distribution import EmpiricalDistribution
from network.finite_hypergraph import FiniteHypergraph, CppFiniteHypergraph
from network.property import BaseNetworkProperty

class FiniteHypergraphTest(unittest.TestCase):
    """Test for the finite clique complex."""

    def setUp(self) -> None:
        cpp_network = CppFiniteHypergraph(
            3, # max_dimension
            [0, 1, 2, 3, 4, 5, 6, 7, 8], # vertices
            [
                [0, 1, 2, 3],
                [1, 2, 4],
                [1, 3, 4],
                [2, 3, 4],
                [0, 5],
                [1, 5],
                [6, 7],
            ], # interactions
            False, # weighted
        )
        self.network = FiniteHypergraph(cpp_network)

    def test_max_dimension(self) -> None:
        self.assertEqual(self.network.max_dimension, 3)

    def test_number_of_simplices(self) -> None:
        self.assertEqual(self.network.num_simplices(0), 9)
        self.assertEqual(self.network.num_simplices(1), 12)
        self.assertEqual(self.network.num_simplices(2), 7)
        self.assertEqual(self.network.num_simplices(3), 1)

    def test_degree_distribution_0_1(self) -> None:
        degree_distribution: EmpiricalDistribution = self.network.calc_base_property(BaseNetworkProperty.vertex_edge_degree_distribution)

        self.assertEqual(
            sorted(degree_distribution.value_sequence.tolist()), [0, 1, 1, 2, 3, 4, 4, 4, 5,])

    def test_degree_distribution_1_2(self) -> None:
        degree_distribution: EmpiricalDistribution = self.network.calc_base_property(BaseNetworkProperty.edge_triangle_degree_distribution)

        self.assertEqual(
            sorted(degree_distribution.value_sequence.tolist()), [0, 0, 0, 2, 2, 2, 2, 2, 2, 3, 3, 3,])

    def test_degree_distribution_2_3(self) -> None:
        degree_distribution: EmpiricalDistribution = self.network.calc_base_property(BaseNetworkProperty.triangle_tetrahedra_degree_distribution)

        self.assertEqual(
            sorted(degree_distribution.value_sequence.tolist()), [0, 0, 0, 1, 1, 1, 1,])

    def test_mean_degree(self) -> None:
        mean_degree = self.network.calc_base_property(BaseNetworkProperty.mean_degree)
        self.assertAlmostEqual(mean_degree, 2.6666, places=3)

    def test_max_degree(self) -> None:
        mean_degree = self.network.calc_base_property(BaseNetworkProperty.max_degree)
        self.assertEqual(mean_degree, 5)

    def test_simplex_dimension_distribution(self) -> None:
        degree_distribution: EmpiricalDistribution = self.network.calc_base_property(BaseNetworkProperty.simplex_dimension_distribution)
        dimension_list = degree_distribution.value_sequence.tolist()
        dimension_list.count(0)

        self.assertEqual(dimension_list.count(0), 9)
        self.assertEqual(dimension_list.count(1), 12)
        self.assertEqual(dimension_list.count(2), 7)
        self.assertEqual(dimension_list.count(3), 1)

    def test_Betti_numbers(self) -> None:
        betti_numbers = self.network.calc_base_property(BaseNetworkProperty.betti_numbers)
        self.assertEqual(betti_numbers, [3, 1, 1])

    def test_betti_numbers_by_component(self) -> None:
        betti_numbers_by_component = self.network.calc_base_property(BaseNetworkProperty.betti_numbers_by_component)
        self.assertEqual(betti_numbers_by_component, [[1, 1, 1], [1, 0, 0], [1, 0, 0]])

    def test_num_of_vertices_by_component(self) -> None:
        num_of_vertices_by_component = self.network.calc_base_property(BaseNetworkProperty.num_of_vertices_by_component)
        self.assertEqual(num_of_vertices_by_component.tolist(), [6, 2, 1])


if __name__ == '__main__':
    unittest.main()
