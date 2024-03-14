#include <random>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "finite_hypergraph_model.h"
#include "finite_network.h"
#include "point.h"
#include "rectangle.h"

namespace py = pybind11;

FiniteHypergraphModel::FiniteHypergraphModel(const py::array_t<double> &parameters_in, const uint32_t seed)
    : Model{seed}, FiniteModel{seed}, HypergraphModel{parameters_in}
{
}

std::tuple<FiniteNetwork, MarkPositionList, MarkPositionList> FiniteHypergraphModel::generate_network() const
{
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(lambda() * torus_size())(random_number_generator_)};
    const auto vertices{create_points(num_of_vertices)};
    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    const auto num_of_interactions{std::poisson_distribution<uint32_t>(lambda_prime() * torus_size())(random_number_generator_)};
    const auto interactions{create_points(num_of_interactions)};
    const auto interaction_mark_position_pairs{convert_to_mark_position_pairs(interactions)};

    const auto connections{generate_connections(vertices, interactions)};
    const auto simplices{create_simplices_from_connections(connections)};
    FiniteNetwork network{max_dimension(), vertex_ids, simplices};

    return std::make_tuple<>(network, vertex_mark_position_pairs, interaction_mark_position_pairs);
}