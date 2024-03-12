#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "finite_network.h"
#include "infinite_network.h"
#include "simplex.h"
#include "typedefs.h"

PYBIND11_MODULE(network, m)
{
    m.doc() = "pybind11 network plugin";

    m.def("filter_simplices",
          &filter_simplices_interface,
          "Filter simplices based on vertices to keep.");

    py::class_<FiniteNetwork>(m, "FiniteNetwork")
        .def(py::init<const Dimension, const VertexList &, const ISimplexList &>())
        .def("calc_betti_numbers", &FiniteNetwork::calc_betti_numbers)
        .def("calc_persistence_pairs", &FiniteNetwork::calc_persistence_pairs)
        .def("create_simplicial_complex", &FiniteNetwork::create_simplicial_complex)
        .def("expand", &FiniteNetwork::expand)
        .def("get_skeleton", &FiniteNetwork::get_skeleton_interface)
        .def("num_vertices", &FiniteNetwork::num_vertices)
        .def("num_simplices", &FiniteNetwork::num_simplices)
        .def("get_simplices", &FiniteNetwork::get_simplices_interface)
        .def("filter", &FiniteNetwork::get_filtered_network)
        .def("reset", &FiniteNetwork::reset)
        .def("calc_simplex_dimension_distribution", &FiniteNetwork::calc_simplex_dimension_distribution)
        .def("calc_facet_dimension_distribution", &FiniteNetwork::calc_facet_dimension_distribution)
        .def("calc_interaction_dimension_distribution", &FiniteNetwork::calc_interaction_dimension_distribution)
        .def("calc_degree_sequence", &FiniteNetwork::calc_degree_sequence)
        .def_property("vertices", &FiniteNetwork::get_vertices_interface, &FiniteNetwork::set_vertices)
        .def_property("max_dimension", &FiniteNetwork::get_max_dimension, &FiniteNetwork::set_max_dimension)
        .def_property("interactions", &FiniteNetwork::get_interactions_interface, &FiniteNetwork::set_interactions)
        .def("calc_facets", &FiniteNetwork::get_facets_interface);

    py::class_<InfiniteNetwork>(m, "InfiniteNetwork")
        .def(py::init<const Dimension, const VertexList &, const ISimplexList &, const VertexId>())
        .def("create_simplicial_complex", &InfiniteNetwork::create_simplicial_complex)
        .def("expand", &InfiniteNetwork::expand)
        .def("get_skeleton", &InfiniteNetwork::get_skeleton_interface)
        .def("num_vertices", &InfiniteNetwork::num_vertices)
        .def("num_simplices", &InfiniteNetwork::num_simplices)
        .def("get_simplices", &InfiniteNetwork::get_simplices_interface)
        .def("filter", &InfiniteNetwork::get_filtered_network)
        .def("reset", &InfiniteNetwork::reset)
        .def("calc_simplex_dimension_distribution", &InfiniteNetwork::calc_simplex_dimension_distribution)
        .def("calc_facet_dimension_distribution", &InfiniteNetwork::calc_facet_dimension_distribution)
        .def("calc_interaction_dimension_distribution", &InfiniteNetwork::calc_interaction_dimension_distribution)
        .def("calc_degree_sequence", &InfiniteNetwork::calc_degree_sequence)
        .def_property("vertices", &InfiniteNetwork::get_vertices_interface, &InfiniteNetwork::set_vertices)
        .def_property("max_dimension", &InfiniteNetwork::get_max_dimension, &InfiniteNetwork::set_max_dimension)
        .def_property("interactions", &InfiniteNetwork::get_interactions_interface, &InfiniteNetwork::set_interactions)
        .def("calc_facets", &InfiniteNetwork::get_facets_interface);
}
