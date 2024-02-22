#include <iostream>
#include <limits>
#include <numeric>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "connections_common.h"
#include "connections_infinite.h"
#include "globals.h"
#include "interface.h"
#include "model_parameters.h"
#include "numpy_cpp_conversion.h"
#include "point.h"
#include "rectangle.h"

namespace py = pybind11;

NetworkListInterface generate_infinite_networks_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t num_of_infinite_networks,
    const uint32_t seed)
{
    random_number_generator.seed(seed);
    auto model_parameters{ModelParameters(model_parameters_input)};
    const auto beta{model_parameters.beta};
    const auto gamma{model_parameters.gamma};
    const auto gamma_prime{model_parameters.gamma_prime};
    const auto interaction_intensity{model_parameters.interaction_intensity};

    std::uniform_real_distribution<Mark> mark_distribution(0., 1.);
    NetworkListInterface result{};

    for (auto network_index{0U}; network_index < num_of_infinite_networks; ++network_index)
    {
        const auto u{mark_distribution(random_number_generator)}; // birth time of the typical node
        // expected number of connecting interactions: 2 * b * u^(-g) / (1 - g')

        std::poisson_distribution<> poisson_distribution_interactions(2 * interaction_intensity * beta * std::pow(u, -gamma) / (1 - gamma_prime));
        const auto num_of_interactions{poisson_distribution_interactions(random_number_generator)};
        const auto interactions{create_interactions_infinite(u, beta, gamma, gamma_prime, interaction_intensity)};
        auto vertices{create_vertices_infinite(interactions, beta, gamma)};
        vertices.push_back(Point(0, u, 0., gamma)); // add the typical node

        auto interaction_rectangles{create_rectangles(interactions, gamma_prime)};
        auto vertex_rectangles{create_rectangles(vertices, gamma)};

        const auto connections{generate_connections(
            interaction_rectangles,
            vertex_rectangles,
            model_parameters.beta,
            std::numeric_limits<Position>::infinity())};
        const auto return_value_1{to_numpy(connections)};

        MarkPositionList interaction_mark_position_pairs{};
        interaction_mark_position_pairs.reserve(num_of_interactions);
        for (const auto &interaction : interactions)
        {
            interaction_mark_position_pairs.push_back(std::make_pair(interaction.mark(), interaction.position()));
        }
        const auto return_value_2{to_numpy(interaction_mark_position_pairs)};

        MarkPositionList vertex_mark_position_pairs{};
        vertex_mark_position_pairs.reserve(vertices.size());
        for (const auto &vertex : vertices)
        {
            vertex_mark_position_pairs.push_back(std::make_pair(vertex.mark(), vertex.position()));
        }
        const auto return_value_3{to_numpy(vertex_mark_position_pairs)};

        result.push_back(std::tuple(return_value_1, return_value_2, return_value_3));
        std::cout << "\rGenerating infinite networks: " << network_index + 1U << "/" << num_of_infinite_networks;
    }
    return result;
}

PointList create_interactions_infinite(
    const Mark u,
    const double beta,
    const double gamma,
    const double gamma_prime,
    const double interaction_intensity)
{
    // expected number of connecting interactions: 2 * b * u^(-g) / (1 - g')
    std::poisson_distribution<uint32_t> poisson_distribution_interactions(2 * interaction_intensity * beta * std::pow(u, -gamma) / (1 - gamma_prime));
    const auto num_of_interactions{poisson_distribution_interactions(random_number_generator)};
    const auto interaction_marks{generate_marks(num_of_interactions)};
    const auto interaction_positions{generate_positions_infinite(
        interaction_marks, beta * std::pow(u, -gamma), gamma_prime)};

    PointList interactions{};
    interactions.reserve(num_of_interactions);
    for (auto index{0U}; index < num_of_interactions; ++index)
    {
        interactions.push_back(Point(index, interaction_marks[index], interaction_positions[index], gamma_prime));
    }

    return interactions;
}

PointList create_vertices_infinite(
    const PointList &interactions,
    const double beta,
    const double gamma)
{
    PointList vertices{};
    auto vertex_id{1U}; // 0 is reserved for the typical vertex
    for (auto interaction_index{0U}; interaction_index < interactions.size(); ++interaction_index)
    {
        const auto vertices_in_neighborhood{
            create_points_in_neighborhood(interactions[interaction_index], beta, gamma)};
        for (const auto &potential_vertex : vertices_in_neighborhood)
        {
            auto should_be_discarded{false};
            for (auto index{0U}; index < interaction_index; ++index)
            {
                if (potential_vertex.connects(interactions[index], beta, std::numeric_limits<float>::infinity()))
                {
                    should_be_discarded = true;
                    break;
                }
            }
            if (!should_be_discarded)
            {
                vertices.push_back(Point(vertex_id, potential_vertex.mark(), potential_vertex.position(), gamma));
                ++vertex_id;
            }
        }
    }
    return vertices;
}

PointList create_points_in_neighborhood(const Point &point, const double beta, const double other_exponent)
{
    PointList points{};
    const auto expected_num_of_points{2 * beta * point.mark_to_gamma() / (1 - other_exponent)};
    const auto num_of_points{std::poisson_distribution<uint32_t>(expected_num_of_points)(random_number_generator)};
    const auto marks{generate_marks(num_of_points)};
    const auto positions{generate_positions_infinite(
        marks, beta * point.mark_to_gamma(), other_exponent)};
    for (auto index{0U}; index < num_of_points; ++index)
    {
        points.push_back(Point(index, marks[index], positions[index], other_exponent));
    }
    return points;
}

PositionList generate_positions_infinite(
    const MarkList &marks,
    const float beta_x_mark_to_gamma,
    const float exponent)
{
    PositionList positions{};
    positions.reserve(marks.size());
    std::uniform_real_distribution<Position> uniform_distribution(-1., 1.);
    for (const auto mark : marks)
    {
        const auto position{uniform_distribution(random_number_generator) * beta_x_mark_to_gamma * std::pow(mark, -exponent)};
        positions.push_back(position);
    }
    return positions;
}
