#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "network.h"
#include "typedefs.h"

PYBIND11_MODULE(network, m)
{
    m.doc() = "pybind11 critical_sections plugin";
    py::class_<Network>(m, "Network")
        .def(py::init<const Dimension>())
        .def("add_simplices", &Network::add_simplices)
        .def("calc_betti_numbers", &Network::calc_betti_numbers)
        .def("calc_persistence_pairs", &Network::calc_persistence_pairs)
        .def("expand", &Network::expand)
        .def("get_skeleton", &Network::get_skeleton)
        .def("num_vertices", &Network::num_vertices)
        .def("num_simplices", &Network::num_simplices)
        .def("keep_only_vertices", &Network::keep_only_vertices)
        .def("reset", &Network::reset)
        .def("calc_simplex_dimension_distribution", &Network::calc_simplex_dimension_distribution)
        .def("calc_degree_sequence", &Network::calc_degree_sequence)
        .def_property("max_dimension", &Network::get_max_dimension, &Network::set_max_dimension)
        .def_property("simplices", &Network::get_simplices, &Network::set_simplices)
        .def_property("facets", &Network::get_facets, &Network::set_facets);
}
