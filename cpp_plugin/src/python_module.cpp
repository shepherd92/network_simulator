#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "finite_adrcm_model.hpp"
#include "finite_clique_complex.hpp"
#include "finite_hypergraph_model.hpp"
#include "finite_hypergraph.hpp"
#include "infinite_adrcm_model.hpp"
#include "infinite_hypergraph_model.hpp"
#include "infinite_hypergraph.hpp"
#include "infinite_clique_complex.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cpp_plugin, m)
{
    m.doc() = "pybind11 model plugin";

    py::class_<FiniteAdrcmModel>(m, "FiniteAdrcmModel")
        .def(py::init<const std::vector<double> &, const uint32_t>())
        .def("generate_network", &FiniteAdrcmModel::generate_network);

    py::class_<InfiniteAdrcmModel>(m, "InfiniteAdrcmModel")
        .def(py::init<const std::vector<double> &, const uint32_t>())
        .def("generate_networks", &InfiniteAdrcmModel::generate_networks);

    py::class_<FiniteHypergraphModel>(m, "FiniteHypergraphModel")
        .def(py::init<const std::vector<double> &, const uint32_t>())
        .def("generate_network", &FiniteHypergraphModel::generate_network);

    py::class_<InfiniteHypergraphModel>(m, "InfiniteHypergraphModel")
        .def(py::init<const std::vector<double> &, const uint32_t>())
        .def("generate_networks", &InfiniteHypergraphModel::generate_networks)
        .def("generate_network", &InfiniteHypergraphModel::generate_network)
        .def("calc_vertex_interaction_degree_sequence_directly",
             &InfiniteHypergraphModel::calc_vertex_interaction_degree_sequence_directly);

    py::class_<Network>(m, "Network")
        .def("get_skeleton", &Network::get_skeleton_interface)
        .def("calc_simplex_dimension_distribution", &Network::calc_simplex_dimension_distribution)
        .def("num_vertices", &Network::num_vertices)
        .def("num_simplices", &Network::num_simplices)
        .def("get_simplices", &Network::get_simplices)
        .def_property_readonly("max_dimension", &Network::get_max_dimension)
        .def_property("vertices", &Network::get_vertices, &Network::set_vertices);

    py::class_<FiniteNetwork, Network>(m, "FiniteNetwork")
        .def("calc_coface_degree_sequence", &FiniteNetwork::calc_coface_degree_sequence)
        .def("calc_betti_numbers", &FiniteNetwork::calc_betti_numbers);

    py::class_<InfiniteNetwork, Network>(m, "InfiniteNetwork")
        .def("calc_coface_degree_sequence", &InfiniteNetwork::calc_coface_degree_sequence);

    py::class_<Hypergraph, Network>(m, "Hypergraph")
        .def("calc_interaction_dimension_distribution", &Hypergraph::calc_interaction_dimension_distribution)
        .def_property("interactions", &Hypergraph::get_interactions, &Hypergraph::set_interactions);

    py::class_<FiniteCliqueComplex, FiniteNetwork, Network>(m, "FiniteCliqueComplex")
        .def(py::init<const Dimension, const PointIdList &, const ConnectionList &>())
        .def("filter", &FiniteCliqueComplex::filter)
        .def_property("edges", &FiniteCliqueComplex::get_edges, &FiniteCliqueComplex::set_edges);

    py::class_<InfiniteCliqueComplex, InfiniteNetwork, Network>(m, "InfiniteCliqueComplex")
        .def(py::init<const Dimension, const PointIdList &, const Mark, const MarkList &>())
        .def("filter", &InfiniteCliqueComplex::filter);

    py::class_<FiniteHypergraph, FiniteNetwork, Hypergraph, Network>(m, "FiniteHypergraph")
        .def(py::init<const Dimension, const PointIdList &, const ISimplexList &, const bool>())
        .def("calc_persistence_intervals", &FiniteHypergraph::calc_persistence_intervals)
        .def("calc_persistence_pairs", &FiniteHypergraph::calc_persistence_pairs)
        .def("calc_simplex_interaction_degree_sequence", &FiniteHypergraph::calc_simplex_interaction_degree_sequence)
        .def("calc_vertex_interaction_degree_distribution", &FiniteHypergraph::calc_vertex_interaction_degree_distribution)
        .def("filter", &FiniteHypergraph::filter);

    py::class_<InfiniteHypergraph, InfiniteNetwork, Hypergraph, Network>(m, "InfiniteHypergraph")
        .def(py::init<const Dimension, const PointIdList &, const ISimplexList &, const Mark, const MarkList &>())
        .def("calc_interaction_dimension_distribution", &InfiniteHypergraph::calc_interaction_dimension_distribution)
        .def("calc_simplex_interaction_degree_sequence", &InfiniteHypergraph::calc_simplex_interaction_degree_sequence)
        .def("calc_vertex_interaction_degree_distribution", &InfiniteHypergraph::calc_vertex_interaction_degree_distribution)
        .def("filter", &InfiniteHypergraph::filter);
}
