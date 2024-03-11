#include <iostream>
#include <numeric>
#include <random>
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

NetworkInterface generate_finite_network_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed)
{
    random_number_generator.seed(seed);
    const auto model_parameters{ModelParameters(model_parameters_input)};
    const auto torus_size{model_parameters.network_size};
    const auto num_of_interactions{std::poisson_distribution<uint32_t>(
        model_parameters.interaction_intensity * model_parameters.network_size)(random_number_generator)};
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(
        model_parameters.network_size)(random_number_generator)};

    const auto interactions{create_points_finite(num_of_interactions, torus_size, model_parameters.gamma_prime)};
    auto interaction_rectangles{create_rectangles(interactions, model_parameters.gamma_prime)};

    const auto vertices{create_points_finite(num_of_vertices, torus_size, model_parameters.gamma)};
    auto vertex_rectangles{create_rectangles(vertices, model_parameters.gamma)};

    const auto connections{generate_connections(
        interaction_rectangles,
        vertex_rectangles,
        model_parameters.beta,
        torus_size)};

    const auto return_value_1{to_numpy(connections)};

    MarkPositionList interaction_mark_position_pairs{};
    interaction_mark_position_pairs.reserve(num_of_interactions);

    std::transform(
        interactions.begin(),
        interactions.end(),
        std::back_inserter(interaction_mark_position_pairs),
        [](const auto &interaction)
        {
            return std::make_pair(interaction.mark(), interaction.position());
        });
    const auto return_value_2{to_numpy(interaction_mark_position_pairs)};

    MarkPositionList vertex_mark_position_pairs{};
    vertex_mark_position_pairs.reserve(num_of_vertices);
    std::transform(
        vertices.begin(),
        vertices.end(),
        std::back_inserter(vertex_mark_position_pairs),
        [](const auto &vertex)
        {
            return std::make_pair(vertex.mark(), vertex.position());
        });
    const auto return_value_3{to_numpy(vertex_mark_position_pairs)};

    return std::make_tuple(std::move(return_value_1), std::move(return_value_2), std::move(return_value_3));
}

PointList create_points_finite(const size_t num_of_nodes, const double torus_size, const float exponent)
{
    const auto marks{generate_marks(num_of_nodes)};
    const auto positions{generate_positions_finite(num_of_nodes, torus_size)};

    std::uniform_real_distribution<> uniform_distribution_u(0., 1.);

    PointList nodes{};
    nodes.reserve(num_of_nodes);
    for (auto index{0U}; index < num_of_nodes; ++index)
    {
        nodes.push_back(Point(index, marks[index], positions[index], exponent));
    }

    return nodes;
}

PositionList generate_positions_finite(const size_t num_of_nodes, const double torus_size)
{
    std::uniform_real_distribution<Position> uniform_distribution(
        -torus_size / 2., +torus_size / 2.);

    PositionList positions(num_of_nodes, 0.0);
    for (auto i{0U}; i < num_of_nodes; ++i)
    {
        positions.at(i) = uniform_distribution(random_number_generator);
    }

    return positions;
}
