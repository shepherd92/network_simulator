#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "connection_generator.h"

PYBIND11_MODULE(adrcm, m)
{
    m.doc() = "pybind11 ADRCM plugin";
    m.def(
        "generate_finite_network_connections",
        &generate_finite_network_connections_interface,
        "Generate connections for the age dependent random simplex model.");
    m.def(
        "generate_infinite_network_connections",
        &generate_infinite_network_connections_interface,
        "Generate connections for the age dependent random simplex model.");
}
