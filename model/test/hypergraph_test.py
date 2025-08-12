#!/usr/bin/env python3
"""Unittest for AdrcmModel."""

import unittest

from model.hypergraph import HypergraphModel
from network.property import BaseNetworkProperty

class HypergraphModelTest(unittest.TestCase):
    def setUp(self) -> None:
        parameters = HypergraphModel.Parameters(
            max_dimension=2,
            network_size=10.,
            alpha=0.5,
            beta=1.,
            gamma=0.5
        )
        self.model = HypergraphModel()
        self.model.parameters = parameters
        self.seed = 1

    def test_generate_finite_network(self) -> None:
        network = self.model.generate_finite_network(self.seed)
        self.assertEqual(network.max_dimension, 2, 'Incorrect max_dimension')
        self.assertEqual(network.num_simplices(0), 13, 'Incorrect number of vertices')
        self.assertEqual(network.num_simplices(1), 78, 'Incorrect number of edges')

    def test_generate_infinite_network_set(self) -> None:
        network_set = self.model.generate_infinite_network_set(10, self.seed)
        self.assertEqual(len(network_set._infinite_networks), 10)
        largest_network = network_set.get_largest_network()
        self.assertEqual(largest_network.max_dimension, 2)


if __name__ == '__main__':
    unittest.main()
