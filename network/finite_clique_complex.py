#!/usr/bin/env python3
"""Finite clique complex network type."""

from typing import Any, Self

from network.finite_network import FiniteNetwork

# pylint: disable-next=no-name-in-module
from cpp_plugin.build.release.cpp_plugin import FiniteCliqueComplex as CppFiniteCliqueComplex


class FiniteCliqueComplex(FiniteNetwork):
    """Base class representing a data set or a simulated finite network."""

    def info(self) -> dict[str, Any]:
        """Return a dict representation based on the network properties."""
        return super().info() | {
            'num_of_edges': self.num_simplices(1),
        }

    def _construct_self_from_cpp(self, cpp_network: CppFiniteCliqueComplex) -> Self:
        """Construct self from the given cpp network."""
        return FiniteCliqueComplex(cpp_network)
