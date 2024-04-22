#!/usr/bin/env python3
"""Levy alpha-stable theoretical distribution test."""

from pathlib import Path
import unittest

import numpy as np

from data_set.data_set import DataSet
from data_set.test_data_set import TestDataSet
from distribution.empirical_distribution import EmpiricalDistribution
from network.property import BaseNetworkProperty


class DataSetTest(unittest.TestCase):
    """Test case for testing the StableDistribution class."""

    def setUp(self):
        """Set up test case."""
        data_set_properties = DataSet.Parameters(
            location=Path(),
            max_dimension=3,
            max_simplex_dimension=3,
            component_index_from_largest=-1,
        )
        self.data_set = TestDataSet(data_set_properties)
        self.data_set.load_data()

    def test_degree_distribution(self):
        """Test if the higher order degree distribution is correct."""
        degree_distribution: EmpiricalDistribution = \
            self.data_set.calc_base_property(BaseNetworkProperty.vertex_edge_degree_distribution)

        np.testing.assert_array_equal(
            degree_distribution.value_sequence,
            np.array([0, 2, 2, 3, 4, 4, 4, 4, 5, ])
        )

    def test_ho_degree_distribution(self):
        """Test if the higher order degree distribution is correct."""
        degree_distribution: EmpiricalDistribution = \
            self.data_set.calc_base_property(BaseNetworkProperty.edge_triangle_degree_distribution)

        np.testing.assert_array_equal(
            degree_distribution.value_sequence,
            np.array([0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, ])
        )


if __name__ == '__main__':

    unittest.main()
