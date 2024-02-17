#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "connections_common.h"
#include "connections_finite.h"
#include "globals.h"
#include "interface.h"
#include "model_parameters.h"
#include "numpy_cpp_conversion.h"
#include "point.h"
#include "rectangle.h"

namespace py = pybind11;

std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>> generate_finite_network_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed)
{
    const auto model_parameters{ModelParameters(model_parameters_input)};

    random_number_generator.seed(seed);

    const auto interactions{create_points_finite(model_parameters.num_of_interactions, model_parameters.torus_size, model_parameters.gamma_prime)};
    auto interaction_rectangles{create_rectangles(
        model_parameters.torus_size,
        model_parameters.num_of_interactions,
        model_parameters.gamma_prime)};
    fill_rectangles(interaction_rectangles, interactions);

    const auto nodes{create_points_finite(model_parameters.num_of_nodes, model_parameters.torus_size, model_parameters.gamma)};
    auto vertex_rectangles{create_rectangles(
        model_parameters.torus_size,
        model_parameters.num_of_nodes,
        model_parameters.gamma)};
    fill_rectangles(vertex_rectangles, nodes);

    const auto connections{generate_connections(interaction_rectangles, vertex_rectangles, true, model_parameters)};

    const auto return_value_1{vector_of_pairs_to_numpy<int>(connections)};

    std::vector<std::pair<float, float>> interaction_mark_position_pairs{};
    interaction_mark_position_pairs.reserve(model_parameters.num_of_interactions);
    for (const auto &interaction : interactions)
    {
        interaction_mark_position_pairs.push_back(std::make_pair(interaction.mark(), interaction.position()));
    }
    const auto return_value_2{vector_of_pairs_to_numpy<float>(interaction_mark_position_pairs)};

    std::vector<std::pair<float, float>> node_mark_position_pairs{};
    node_mark_position_pairs.reserve(model_parameters.num_of_nodes);
    for (const auto &node : nodes)
    {
        node_mark_position_pairs.push_back(std::make_pair(node.mark(), node.position()));
    }
    const auto return_value_3{vector_of_pairs_to_numpy<float>(node_mark_position_pairs)};

    return std::make_tuple(std::move(return_value_1), std::move(return_value_2), std::move(return_value_3));
}

std::vector<Point> create_points_finite(const size_t num_of_nodes, const double torus_size, const float exponent)
{
    const auto marks{generate_marks(num_of_nodes)};
    const auto positions{generate_positions_finite(num_of_nodes, torus_size)};

    std::uniform_real_distribution<> uniform_distribution_u(0., 1.);

    std::vector<Point> nodes{};
    nodes.reserve(num_of_nodes);
    for (auto index{0U}; index < num_of_nodes; ++index)
    {
        nodes.push_back(Point(index, marks[index], positions[index], exponent));
    }

    return nodes;
}

std::vector<float> generate_positions_finite(
    const size_t num_of_nodes,
    const double torus_size)
{
    std::uniform_real_distribution<float> uniform_distribution(
        -torus_size / 2., +torus_size / 2.);

    std::vector<float> positions(num_of_nodes, 0.0);
    for (auto i{0U}; i < num_of_nodes; ++i)
    {
        positions.at(i) = uniform_distribution(random_number_generator);
    }

    return positions;
}
