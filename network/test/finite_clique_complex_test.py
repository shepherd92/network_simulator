#!/usr/bin/env python3
"""Unittest for AdrcmModel."""

import unittest

from distribution.empirical_distribution import EmpiricalDistribution
from network.finite_clique_complex import FiniteCliqueComplex, CppFiniteCliqueComplex
from network.property import BaseNetworkProperty

class FiniteCliqueComplexTest(unittest.TestCase):
    """Test for the finite clique complex."""

    def setUp(self) -> None:
        cpp_network = CppFiniteCliqueComplex(
            3, # max_dimension
            [0, 1, 2, 3, 4, 5, 6, 7], # vertices
            [
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 2),
                (1, 3),
                (2, 3),
                (2, 4),
                (3, 4),
                (5, 6),
            ], # edges
        )
        self.network = FiniteCliqueComplex(cpp_network)

    def test_number_of_simplices(self) -> None:
        self.assertEqual(self.network.num_simplices(2), 5)

    def test_degree_distribution_0_1(self) -> None:
        degree_distribution: EmpiricalDistribution = self.network.calc_base_property(BaseNetworkProperty.vertex_edge_degree_distribution)

        self.assertEqual(
            sorted(degree_distribution.value_sequence), [0, 1, 1, 2, 3, 3, 4, 4])

    def test_degree_distribution_1_2(self) -> None:
        degree_distribution: EmpiricalDistribution = self.network.calc_base_property(BaseNetworkProperty.edge_triangle_degree_distribution)

        self.assertEqual(
            sorted(degree_distribution.value_sequence), [0, 1, 1, 2, 2, 2, 2, 2, 3])

    def test_simplex_dimension_distribution(self) -> None:
        degree_distribution: EmpiricalDistribution = self.network.calc_base_property(BaseNetworkProperty.simplex_dimension_distribution)
        dimension_list = degree_distribution.value_sequence.tolist()
        dimension_list.count(0)

        self.assertEqual(dimension_list.count(0), 8)
        self.assertEqual(dimension_list.count(1), 9)
        self.assertEqual(dimension_list.count(2), 5)
        self.assertEqual(dimension_list.count(3), 1)

    def test_Betti_numbers(self) -> None:
        betti_numbers = self.network.calc_base_property(BaseNetworkProperty.betti_numbers)
        self.assertEqual(betti_numbers, [3, 0, 0])


if __name__ == '__main__':
    unittest.main()
