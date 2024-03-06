#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "skeleton_blocker_finite_network.h"
#include "skeleton_blocker_infinite_network.h"
#include "simplex_tree_finite_network.h"
#include "simplex_tree_infinite_network.h"
#include "typedefs.h"

PYBIND11_MODULE(network, m)
{
    m.doc() = "pybind11 network plugin";

    py::class_<SimplexTreeFiniteNetwork>(m, "SimplexTreeFiniteNetwork")
        .def(py::init<const Dimension>())
        .def("add_simplices", &Network::add_simplices_interface)
        .def("calc_betti_numbers", &SimplexTreeFiniteNetwork::calc_betti_numbers)
        .def("calc_persistence_pairs", &SimplexTreeFiniteNetwork::calc_persistence_pairs)
        .def("expand", &SimplexTreeNetwork::expand)
        .def("get_skeleton", &Network::get_skeleton_interface)
        .def("num_vertices", &Network::num_vertices)
        .def("num_simplices", &Network::num_simplices)
        .def("keep_only_vertices", &Network::keep_only_vertices)
        .def("reset", &Network::reset)
        .def("calc_simplex_dimension_distribution", &Network::calc_simplex_dimension_distribution)
        .def("calc_facet_dimension_distribution", &Network::calc_facet_dimension_distribution)
        .def("calc_interaction_dimension_distribution", &Network::calc_interaction_dimension_distribution)
        .def("calc_degree_sequence", &Network::calc_degree_sequence)
        .def_property("max_dimension", &SimplexTreeNetwork::get_max_dimension, &SimplexTreeNetwork::set_max_dimension)
        .def_property("vertices", &Network::get_vertices_interface, &Network::set_vertices)
        .def_property("interactions", &Network::get_interactions_interface, &Network::set_interactions)
        .def_property("facets", &Network::get_facets_interface, &Network::set_facets);

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
