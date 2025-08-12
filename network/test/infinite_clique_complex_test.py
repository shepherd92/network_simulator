#!/usr/bin/env python3
"""Unittest for InfiniteCliqueComplex."""

import unittest

from network.infinite_clique_complex import InfiniteCliqueComplex, CppInfiniteCliqueComplex
from network.property import BaseNetworkProperty

class InfiniteCliqueComplexTest(unittest.TestCase):
    """Test for the finite clique complex."""

    def setUp(self) -> None:
        cpp_network = CppInfiniteCliqueComplex(
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
            0.45, # typical_mark
            [0.7, 0.5, 0.5, 0.6, 0.8, 0.3, 0.9, 0.8], # mark_list
        )
        self.network = InfiniteCliqueComplex(cpp_network)

    def test_number_of_simplices(self) -> None:
        # typical vertex is not included in the simplices
        self.assertEqual(self.network.num_simplices(0), 0)
        self.assertEqual(self.network.num_simplices(1), 7)
        self.assertEqual(self.network.num_simplices(2), 8)
        self.assertEqual(self.network.num_simplices(3), 5)

    def test_degree_distribution_0_1(self) -> None:

        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.vertex_edge_degree_distribution)
        self.assertEqual(sorted(degrees), [8])

    def test_degree_distribution_1_2(self) -> None:

        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.edge_triangle_degree_distribution)
        self.assertEqual(sorted(degrees), [0, 1, 2, 3, 3, 4, 4])

    def test_degree_distribution_2_3(self) -> None:

        degrees: list[int] = self.network.calc_base_property(BaseNetworkProperty.triangle_tetrahedra_degree_distribution)
        self.assertEqual(sorted(degrees), [1, 1, 2, 2, 2, 2, 2, 3])


if __name__ == '__main__':
    unittest.main()
