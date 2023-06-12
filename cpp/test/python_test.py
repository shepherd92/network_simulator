#!/usr/bin/env python3
"""C++ sections test."""

import unittest

from cpp.network import Network


class CppModulesTest(unittest.TestCase):
    """Test case for testing the StableDistribution class."""

    def setUp(self):
        """Set up test case."""
        self.network = Network()

    def test_get_facets(self):
        """Test if the fitting method gives a reasonably good fit."""
        facets = self.network.get_facets()
        return


if __name__ == '__main__':

    unittest.main()
