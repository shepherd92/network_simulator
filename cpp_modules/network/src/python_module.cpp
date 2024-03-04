#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "finite_network.h"
#include "infinite_network.h"
#include "typedefs.h"

PYBIND11_MODULE(network, m)
{
    m.doc() = "pybind11 network plugin";

    py::class_<FiniteNetwork>(m, "FiniteNetwork")
        .def(py::init<const Dimension>())
        .def("add_simplices", &FiniteNetwork::add_simplices_interface)
        .def("calc_betti_numbers", &FiniteNetwork::calc_betti_numbers)
        .def("calc_persistence_pairs", &FiniteNetwork::calc_persistence_pairs)
        .def("expand", &FiniteNetwork::expand)
        .def("get_skeleton", &FiniteNetwork::get_skeleton_interface)
        .def("num_vertices", &FiniteNetwork::num_vertices)
        .def("num_simplices", &FiniteNetwork::num_simplices)
        .def("keep_only_vertices", &FiniteNetwork::keep_only_vertices)
        .def("reset", &FiniteNetwork::reset)
        .def("calc_simplex_dimension_distribution", &FiniteNetwork::calc_simplex_dimension_distribution)
        .def("calc_facet_dimension_distribution", &FiniteNetwork::calc_facet_dimension_distribution)
        .def("calc_interaction_dimension_distribution", &FiniteNetwork::calc_interaction_dimension_distribution)
        .def("calc_degree_sequence", &FiniteNetwork::calc_degree_sequence)
        .def_property("max_dimension", &FiniteNetwork::get_max_dimension, &FiniteNetwork::set_max_dimension)
        .def_property("vertices", &FiniteNetwork::get_vertices_interface, &FiniteNetwork::set_vertices)
        .def_property("simplices", &FiniteNetwork::get_simplices_interface, &FiniteNetwork::set_simplices)
        .def_property("interactions", &FiniteNetwork::get_interactions_interface, &FiniteNetwork::set_interactions)
        .def_property("facets", &FiniteNetwork::get_facets_interface, &FiniteNetwork::set_facets);

    py::class_<InfiniteNetwork>(m, "InfiniteNetwork")
        .def(py::init<const Dimension, const VertexId>())
        .def("add_vertices", &InfiniteNetwork::add_vertices)
        .def("add_simplices", &InfiniteNetwork::add_simplices)
        .def("calc_betti_numbers", &InfiniteNetwork::calc_betti_numbers)
        .def("calc_persistence_pairs", &InfiniteNetwork::calc_persistence_pairs)
        .def("expand", &InfiniteNetwork::expand)
        .def("get_skeleton", &InfiniteNetwork::get_skeleton_interface)
        .def("num_vertices", &InfiniteNetwork::num_vertices)
        .def("num_simplices", &InfiniteNetwork::num_simplices)
        .def("keep_only_vertices", &InfiniteNetwork::keep_only_vertices)
        .def("reset", &InfiniteNetwork::reset)
        .def("calc_simplex_dimension_distribution", &InfiniteNetwork::calc_simplex_dimension_distribution)
        .def("calc_facet_dimension_distribution", &InfiniteNetwork::calc_facet_dimension_distribution)
        .def("calc_interaction_dimension_distribution", &InfiniteNetwork::calc_interaction_dimension_distribution)
        .def("calc_degree_sequence", &InfiniteNetwork::calc_degree_sequence)
        .def_property("max_dimension", &InfiniteNetwork::get_max_dimension, &InfiniteNetwork::set_max_dimension)
        .def_property("vertices", &InfiniteNetwork::get_vertices_interface, &InfiniteNetwork::set_vertices)
        .def_property("simplices", &InfiniteNetwork::get_simplices_interface, &InfiniteNetwork::set_simplices)
        .def_property("interactions", &InfiniteNetwork::get_interactions_interface, &InfiniteNetwork::set_interactions)
        .def_property("facets", &InfiniteNetwork::get_facets_interface, &InfiniteNetwork::set_facets);
}
