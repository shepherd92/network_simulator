#!/usr/bin/env python3
"""C++ sections test."""

import unittest
import numpy as np
from numpy.testing import assert_array_equal

from cpp.build.cpp_library import Network


class CppModulesTest(unittest.TestCase):
    """Test case for testing the StableDistribution class."""

    def setUp(self):
        """Set up test case."""
        self.network = Network(3)

    def test_get_simplices_by_dimension(self) -> None:
        """Test if the fitting method gives a reasonably good fit."""
        self.network.add_simplices([[1], [0, 1, 2], [0, 1], [3, 5, 4], [3, 4, 5, 6]])
        simplices = self.network.get_simplices_by_dimension(2)
        assert_array_equal(simplices, np.array([
            [0, 1, 2], [3, 4, 5], [3, 4, 6], [3, 5, 6], [4, 5, 6]
        ], dtype=np.int32))

    def test_facets(self) -> None:
        """Test if the fitting method gives a reasonably good fit."""
        self.network.add_simplices([[1], [0, 1, 2], [0, 1], [3, 5, 4], [3, 4, 5, 6]])
        facets = self.network.facets()
        assert facets == [[0, 1, 2], [3, 4, 5, 6]]


if __name__ == '__main__':

    unittest.main()
