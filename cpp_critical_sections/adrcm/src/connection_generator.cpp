#include <iostream>
#include <random>
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include "connection_generator.h"
#include "model_parameters.h"
#include "numpy_cpp_conversion.h"
#include "point.h"

namespace py = pybind11;

bool is_close(const double first, const double second);
std::vector<Point> create_nodes(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input);
std::vector<std::pair<int, int>> generate_network_connections_default(
    const std::vector<Point> &nodes,
    const ModelParameters &model_parameters);
double profile_function(const double argument, const double alpha);

py::array_t<int> generate_finite_network_connections_default_interface(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input,
    const py::array_t<double> &model_parameters_input)
{
    const auto nodes{create_nodes(birth_times_input, positions_input)};
    const auto model_parameters{ModelParameters(model_parameters_input)};

    const auto connections{generate_network_connections_default(nodes, model_parameters)};
    return vector_of_pairs_to_numpy<int>(connections);
}

std::vector<py::array_t<int>> generate_infinite_network_connections_default_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed)
{
    auto model_parameters{ModelParameters(model_parameters_input)};
    std::vector<py::array_t<int>> result{};
    std::mt19937 random_number_generator{seed};
    std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
    const auto b{model_parameters.beta};
    const auto g{model_parameters.gamma};

    for (auto network_index{0U}; network_index < model_parameters.num_of_nodes; ++network_index)
    {
        const auto u{uniform_distribution(random_number_generator)}; // birth time of oldest node
        model_parameters.torus_size = b / u;
        const auto N{static_cast<uint32_t>(b * (1. / u - 1.))}; // total area of the rectangle

        if (N > 1e5)
        {
            // std::cout << "Warning: " << N << " nodes should be simulated for infinite networks" << std::endl;
            continue;
        }

        // generate nodes
        std::vector<Point> nodes;
        nodes.reserve(N + 1U);
        nodes.push_back(Point(u, 0.));

        std::uniform_real_distribution<> birth_time_distribution{u, 1.0};
        std::uniform_real_distribution<> position_distribution{-b / u, +b / u};
        for (auto node_index{0U}; node_index < N; ++node_index)
        {
            nodes.push_back(Point(
                birth_time_distribution(random_number_generator),
                position_distribution(random_number_generator)));
        }

        // vertices closer to the origin than 1/2 * beta * u^(-gamma) * v^(gamma - 1) are connected to  the origin
        std::vector<Point> nodes_connected_to_o;
        nodes_connected_to_o.push_back(nodes[0]);
        for (const auto &node : nodes)
        {
            if (node.position() < 0.5 * b * std::pow(u, -g) * std::pow(node.birth_time(), g - 1.))
            {
                nodes_connected_to_o.push_back(node);
            }
        }

        const auto connections{generate_network_connections_default(nodes_connected_to_o, model_parameters)};
        result.push_back(vector_of_pairs_to_numpy<int>(connections));
    }

    return result;
}

std::vector<Point> create_nodes(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input)
{
    const auto birth_times{numpy_to_vector_1d<double>(birth_times_input)};
    const auto positions{numpy_to_vector_1d<double>(positions_input)};
    const auto num_of_nodes{birth_times.size()};
    std::vector<Point> nodes;
    nodes.reserve(num_of_nodes);

    for (auto index{0U}; index < num_of_nodes; ++index)
    {
        nodes.emplace_back(Point(birth_times[index], positions[index]));
    }

    return nodes;
}

std::vector<std::pair<int, int>> generate_network_connections_default(
    const std::vector<Point> &nodes,
    const ModelParameters &model_parameters)
{
    const auto num_of_nodes{nodes.size()};
    std::vector<std::pair<int, int>> connections{};
    for (auto source_node_id{0U}; source_node_id < num_of_nodes; ++source_node_id)
    {
        const auto &source_node{nodes[source_node_id]};

        for (auto target_node_id{0U}; target_node_id < source_node_id; ++target_node_id)
        {
            const auto &target_node{nodes[target_node_id]};
            const auto birth_time_ratio{source_node.birth_time() / target_node.birth_time()};
            const auto distance{source_node.torus_distance(target_node, model_parameters.torus_size)};
            const auto profile_function_argument{
                distance * source_node.birth_time() /
                (model_parameters.beta * pow(birth_time_ratio, model_parameters.gamma))};
            if (profile_function_argument < 0.5)
            {
                connections.push_back(std::pair(source_node_id, target_node_id));
            }
        }
    }
    return connections;
}