#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "connection_generator.h"

PYBIND11_MODULE(hypergraph, m)
{
    m.doc() = "pybind11 hypergraph plugin";
    m.def(
        "generate_finite_network_connections",
        &generate_finite_network_connections_interface,
        "Generate connections for the hypergraph model.");
    m.def(
        "generate_infinite_network_connections",
        &generate_infinite_network_connections_interface,
        "Generate connections for the hypergraph model.");
}
