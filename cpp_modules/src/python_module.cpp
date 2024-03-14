#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "finite_adrcm_model.h"
#include "finite_hypergraph_model.h"
#include "finite_network.h"
#include "infinite_adrcm_model.h"
#include "infinite_hypergraph_model.h"
#include "infinite_network.h"

PYBIND11_MODULE(cpp_plugin, m)
{
    m.doc() = "pybind11 model plugin";

    m.def("filter_simplices",
          &filter_simplices_interface,
          "Filter simplices based on vertices to keep.");

    py::class_<FiniteAdrcmModel>(m, "FiniteAdrcmModel")
        .def(py::init<const py::array_t<double> &, const uint32_t>())
        .def("generate_network", &FiniteAdrcmModel::generate_network);

    py::class_<InfiniteAdrcmModel>(m, "InfiniteAdrcmModel")
        .def(py::init<const py::array_t<double> &, const uint32_t>())
        .def("generate_networks", &InfiniteAdrcmModel::generate_networks);

    py::class_<FiniteHypergraphModel>(m, "FiniteHypergraphModel")
        .def(py::init<const py::array_t<double> &, const uint32_t>())
        .def("generate_network", &FiniteHypergraphModel::generate_network);

    py::class_<InfiniteHypergraphModel>(m, "InfiniteHypergraphModel")
        .def(py::init<const py::array_t<double> &, const uint32_t>())
        .def("generate_networks", &InfiniteHypergraphModel::generate_networks);

    py::class_<FiniteNetwork>(m, "FiniteNetwork")
        .def(py::init<const Dimension, const PointIdList &, const ISimplexList &>())
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
        .def(py::init<const Dimension, const PointIdList &, const ISimplexList &, const PointId>())
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
