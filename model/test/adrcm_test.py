#!/usr/bin/env python3
"""Unittest for AdrcmModel."""

import unittest
from model.adrcm import AdrcmModel

class AdrcmModelTest(unittest.TestCase):
    def setUp(self) -> None:
        parameters = AdrcmModel.Parameters(
            max_dimension=2,
            network_size=10.,
            alpha=0.5,
            beta=1.,
            gamma=0.5
        )
        self.model = AdrcmModel()
        self.model.parameters = parameters

    def test_generate_finite_network(self) -> None:
        seed = 1
        network = self.model.generate_finite_network(seed)
        self.assertEqual(network.max_dimension, 2, 'Incorrect max_dimension')
        self.assertEqual(network.num_simplices(0), 13, 'Incorrect number of vertices')
        self.assertEqual(network.num_simplices(1), 78, 'Incorrect number of edges')

if __name__ == '__main__':
    unittest.main()
