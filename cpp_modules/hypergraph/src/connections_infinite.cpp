#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>
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

std::vector<std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>>> generate_infinite_networks_interface(
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

    std::uniform_real_distribution<float> uniform_distribution(0., 1.);
    std::vector<std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>>> result{};

    for (auto network_index{0U}; network_index < num_of_infinite_networks; ++network_index)
    {
        const auto u{uniform_distribution(random_number_generator)}; // birth time of the typical node
        // expected number of connecting interactions: 2 * b * u^(-g) / (1 - g')

        std::poisson_distribution<> poisson_distribution_interactions(2 * interaction_intensity * beta * std::pow(u, -gamma) / (1 - gamma_prime));
        const auto num_of_interactions{poisson_distribution_interactions(random_number_generator)};
        const auto interactions{create_interactions_infinite(u, beta, gamma, gamma_prime, interaction_intensity)};
        auto vertices{create_vertices_infinite(interactions, beta, gamma, gamma_prime)};
        vertices.push_back(Point(0, u, 0., gamma)); // add the typical node

        const auto network_size{2. * determine_infinite_network_size(interactions, vertices)};
        auto interaction_rectangles{create_rectangles(network_size, interactions.size(), gamma_prime)};
        fill_rectangles(interaction_rectangles, interactions);
        auto vertex_rectangles{create_rectangles(network_size, vertices.size(), gamma)};
        fill_rectangles(vertex_rectangles, vertices);

        const auto connections{generate_connections(interaction_rectangles, vertex_rectangles, false, model_parameters)};
        const auto return_value_1{vector_of_pairs_to_numpy<int>(connections)};

        std::vector<std::pair<float, float>> interaction_mark_position_pairs{};
        interaction_mark_position_pairs.reserve(num_of_interactions);
        for (const auto &interaction : interactions)
        {
            interaction_mark_position_pairs.push_back(std::make_pair(interaction.mark(), interaction.position()));
        }
        const auto return_value_2{vector_of_pairs_to_numpy<float>(interaction_mark_position_pairs)};

        std::vector<std::pair<float, float>> vertex_mark_position_pairs{};
        vertex_mark_position_pairs.reserve(model_parameters.num_of_nodes);
        for (const auto &vertex : vertices)
        {
            vertex_mark_position_pairs.push_back(std::make_pair(vertex.mark(), vertex.position()));
        }
        const auto return_value_3{vector_of_pairs_to_numpy<float>(vertex_mark_position_pairs)};

        result.push_back(std::tuple(return_value_1, return_value_2, return_value_3));
    }
    return result;
}

double determine_infinite_network_size(
    const std::vector<Point> &interactions,
    const std::vector<Point> &vertices)
{
    float max_position{0.};
    for (const auto &interaction : interactions)
    {
        max_position = std::max(max_position, fabs(interaction.position()));
    }
    for (const auto &vertex : vertices)
    {
        max_position = std::max(max_position, fabs(vertex.position()));
    }
    return max_position;
}

std::vector<Point> create_interactions_infinite(
    const double u,
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

    std::vector<Point> interactions{};
    interactions.reserve(num_of_interactions);
    for (auto index{0U}; index < num_of_interactions; ++index)
    {
        interactions.push_back(Point(index, interaction_marks[index], interaction_positions[index], gamma_prime));
    }

    return interactions;
}

std::vector<Point> create_vertices_infinite(
    const std::vector<Point> &interactions,
    const double beta,
    const double gamma,
    const double gamma_prime)
{
    std::vector<Point> vertices{};
    auto vertex_id{1U}; // 0 is reserved for the typical vertex
    for (auto interaction_index{0U}; interaction_index < interactions.size(); ++interaction_index)
    {
        const auto &interaction{interactions[interaction_index]};
        const auto v{interaction.mark()};
        const auto expected_num_of_vertices{2 * beta * std::pow(v, -gamma_prime) / (1 - gamma)};
        const auto num_of_vertices{std::poisson_distribution<>(expected_num_of_vertices)(random_number_generator)};
        const auto vertex_marks{generate_marks(num_of_vertices)};
        const auto vertex_positions{generate_positions_infinite(
            vertex_marks, beta * std::pow(v, -gamma_prime), gamma)};

        for (const auto &vertex : vertices)
        {
            auto should_be_discarded{false};
            for (auto index{0U}; index < interaction_index; ++index)
            {
                if (vertex.connects(interactions[index], 0., beta, false))
                {
                    should_be_discarded = true;
                    break;
                }
            }
            if (!should_be_discarded)
            {
                vertices.push_back(Point(vertex_id, vertex.mark(), vertex.position(), gamma));
                ++vertex_id;
            }
        }
    }
    return vertices;
}

std::vector<float> generate_positions_infinite(
    const std::vector<float> &marks,
    const float beta_x_mark_to_gamma,
    const float exponent)
{
    std::vector<float> positions{};
    positions.reserve(marks.size());
    std::uniform_real_distribution<> uniform_distribution(-1., 1.);
    for (const auto mark : marks)
    {
        const auto position{uniform_distribution(random_number_generator) * beta_x_mark_to_gamma * std::pow(mark, -exponent)};
        positions.push_back(position);
    }
    return positions;
}
