#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "interface.h"

PYBIND11_MODULE(hypergraph, m)
{
    m.doc() = "pybind11 hypergraph plugin";
    m.def(
        "generate_finite_network_cpp",
        &generate_finite_network_interface,
        "Generate connections for the finite hypergraph model.");
    m.def(
        "generate_infinite_networks_cpp",
        &generate_finite_network_interface,
        "Generate connections for the infinite hypergraph model.");
}
