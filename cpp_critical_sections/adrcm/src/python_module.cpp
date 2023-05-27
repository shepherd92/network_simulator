#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "connection_generator.h"

PYBIND11_MODULE(adrcm, m)
{
    m.doc() = "pybind11 ADRCM plugin";
    m.def("generate_connections_default", &generate_connections_default, "Generate connections for the age dependent random simplex model.");
}
