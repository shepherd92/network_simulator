#include <iostream>
#include <numeric>
#include <random>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "connection_generator.h"
#include "globals.h"
#include "model_parameters.h"
#include "numpy_cpp_conversion.h"
#include "point.h"
#include "typedefs.h"

namespace py = pybind11;

bool is_close(const double first, const double second);

PointList create_nodes(const ModelParameters &model_parameters);

std::vector<double> generate_birth_times(const uint32_t num_nodes);
PositionList generate_positions(const ModelParameters &model_parameters);

std::vector<std::pair<int, int>> generate_network_connections(
    const PointList &nodes,
    const bool is_finite,
    const ModelParameters &model_parameters);

py::array_t<int> generate_finite_network_connections_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed)
{
    random_number_generator.seed(seed);
    const auto model_parameters{ModelParameters(model_parameters_input)};
    const auto nodes{create_nodes(model_parameters)};
    const auto connections{generate_network_connections(nodes, true, model_parameters)};

    return to_numpy(connections);
}

std::vector<py::array_t<int>> generate_infinite_network_connections_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t num_of_infinite_networks,
    const uint32_t seed)
{
    random_number_generator.seed(seed);
    auto model_parameters{ModelParameters(model_parameters_input)};
    std::vector<py::array_t<int>> result{};
    const auto a{model_parameters.alpha};
    const auto b{model_parameters.beta};
    const auto g{model_parameters.gamma};
    const auto d{model_parameters.torus_dimension};

    // expected number of incoming connections is: b / g * (u^(-g) - 1)
    std::uniform_real_distribution<> uniform_distribution(0., 1.);
    std::normal_distribution<> normal_distribution(0., 1.);

    for (auto network_index{0U}; network_index < num_of_infinite_networks; ++network_index)
    {
        const auto u{uniform_distribution(random_number_generator)}; // birth time of the typical node

        // generate nodes
        PointList nodes;
        nodes.push_back(Point(u, Position(d, 0.)));

        // generate older nodes which (u, 0) connects to
        // w is a random variable that is later transformed to be the birth time
        // there is a x0.5 due to parameter alpha, but there is a x2 as well due to +-
        const auto w_intensity_older_nodes{a * b / (1. - g) * std::pow(u, g - 1.)};
        std::exponential_distribution<> w_interarrival_time_distribution_older_nodes(w_intensity_older_nodes);
        double w_older{w_interarrival_time_distribution_older_nodes(random_number_generator)};
        while (w_older < std::pow(u, 1. - g))
        {
            // v: birth time of the neighbor of (u, 0)
            const auto v{std::pow(w_older, 1. / (1. - g))};
            const auto radius{std::pow(a * b * std::pow(v, -g) * std::pow(u, g - 1.), 1. / d)};

            // generate point uniformly in the d dimensional sphere (Box-Mueller transform)
            Position position(d);
            auto norm{0.};
            for (auto current_dimension{0U}; current_dimension < d; ++current_dimension)
            {
                const auto gaussian{normal_distribution(random_number_generator)};
                position.push_back(gaussian);
                norm += std::pow(gaussian, 2.);
            }
            norm = std::pow(norm, 0.5);

            // normalize the vector then multiply with U^(1/d) * radius
            for (auto &coordinate : position)
            {
                const auto random_number{uniform_distribution(random_number_generator)};
                coordinate = radius * std::pow(random_number, 1. / d) * coordinate / norm;
            }

            // create point and save it
            nodes.push_back(Point{v, position});
            // increment w to arrive to the next birth time (before transformation)
            w_older += w_interarrival_time_distribution_older_nodes(random_number_generator);
        }

        // generate younger nodes which connect to (u, 0)
        // w is a random variable that is later transformed to be the birth time
        const auto w_intensity_younger_nodes{a * b / g * std::pow(u, -g)};
        std::exponential_distribution<> w_interarrival_time_distribution_younger_nodes{w_intensity_younger_nodes};
        double w_younger{std::pow(u, g) + w_interarrival_time_distribution_younger_nodes(random_number_generator)};
        while (w_younger < 1.)
        {
            // v: birth time of the neighbor of (u, 0)
            const auto v{std::pow(w_younger, 1. / g)};

            const auto radius{std::pow(a * b * std::pow(u, -g) * std::pow(v, g - 1.), 1. / d)};

            // generate point uniformly in the d dimensional sphere (Box-Mueller transform)
            Position position(d);
            auto norm{0.};
            for (auto current_dimension{0U}; current_dimension < d; ++current_dimension)
            {
                const auto gaussian{normal_distribution(random_number_generator)};
                position.push_back(gaussian);
                norm += std::pow(gaussian, 2.);
            }
            norm = std::pow(norm, 0.5);

            // normalize the vector then multiply with U^(1/d) * radius
            for (auto &coordinate : position)
            {
                const auto random_number{uniform_distribution(random_number_generator)};
                coordinate = radius * std::pow(random_number, 1. / d) * coordinate / norm;
            }
            // create point and save it
            nodes.push_back(Point{v, position});
            // increment w to arrive to the next birth time (before transformation)
            w_younger += w_interarrival_time_distribution_younger_nodes(random_number_generator);
        }

        const auto connections{generate_network_connections(nodes, false, model_parameters)};
        result.push_back(to_numpy(connections));
    }
    return result;
}

PointList create_nodes(const ModelParameters &model_parameters)
{
    std::vector<uint32_t> node_ids{model_parameters.num_of_nodes};
    std::iota(node_ids.begin(), node_ids.end(), 0U);
    const auto birth_times{generate_birth_times(model_parameters.num_of_nodes)};
    const auto positions{generate_positions(model_parameters)};

    std::uniform_real_distribution<> uniform_distribution_u(0., 1.);

    PointList nodes{};
    nodes.reserve(model_parameters.num_of_nodes);
    for (auto index{0U}; index < model_parameters.num_of_nodes; ++index)
    {
        nodes.push_back(Point(birth_times[index], positions[index]));
    }

    return nodes;
}

std::vector<double> generate_birth_times(const uint32_t num_nodes)
{
    std::uniform_real_distribution<> uniform_distribution(0., 1.);
    std::vector<double> birth_times{};
    birth_times.reserve(num_nodes);

    for (auto i{0U}; i < num_nodes; ++i)
    {
        birth_times.push_back(uniform_distribution(random_number_generator));
    }

    std::sort(birth_times.begin(), birth_times.end());

    return birth_times;
}

std::vector<std::vector<double>> generate_positions(const ModelParameters &model_parameters)
{
    std::uniform_real_distribution<double> uniform_distribution(
        -model_parameters.torus_size_in_1_dimension / 2.,
        +model_parameters.torus_size_in_1_dimension / 2.);

    std::vector<Position> positions{};
    positions.reserve(model_parameters.num_of_nodes);
    for (auto i{0U}; i < model_parameters.num_of_nodes; ++i)
    {
        Position coordinates{};
        coordinates.reserve(model_parameters.torus_dimension);
        for (auto dimension{0U}; dimension < model_parameters.torus_dimension; ++dimension)
        {
            coordinates.push_back(uniform_distribution(random_number_generator));
        }
        positions.push_back(coordinates);
    }

    return positions;
}

std::vector<std::pair<int, int>> generate_network_connections(
    const PointList &nodes,
    const bool is_finite,
    const ModelParameters &model_parameters)
{
    // create aliases
    const auto a{model_parameters.alpha};
    const auto is_default_a{is_close(a, 0.5)};
    const auto b{model_parameters.beta};
    const auto g{model_parameters.gamma};
    const auto d{model_parameters.torus_dimension};
    const auto num_of_nodes{nodes.size()};
    const auto torus_size{is_finite ? model_parameters.torus_size_in_1_dimension : std::numeric_limits<double>::infinity()};

    std::uniform_real_distribution<> uniform_distribution(0., 1.);
    std::vector<std::pair<int, int>> connections{};
    auto counter{0U};
    for (auto target_node_id{0U}; target_node_id < num_of_nodes; ++target_node_id)
    {
        const auto &target_node{nodes[target_node_id]};
        const auto s{target_node.birth_time()};

        const auto size_of_neighborhood{std::pow((a * b / s), 1. / d)};

        for (auto source_node_id{target_node_id + 1U}; source_node_id < num_of_nodes; ++source_node_id)
        {
            const auto &source_node{nodes[source_node_id]};
            const auto distance{source_node.distance(target_node, torus_size)};
            if (distance > size_of_neighborhood)
            {
                continue;
            }

            const auto t{source_node.birth_time()};
            const auto max_distance_of_connection{std::pow((a * b / t) * std::pow(t / s, g), 1. / d)};

            if (distance < max_distance_of_connection)
            {
                if (is_default_a)
                {
                    connections.push_back(std::pair(source_node_id, target_node_id));
                }
                else
                {
                    const auto random_number{uniform_distribution(random_number_generator)};
                    if (random_number < 1 / (2. * a))
                    {
                        connections.push_back(std::pair(source_node_id, target_node_id));
                    }
                }
            }
        }
        if (++counter % 10000 == 0 && num_of_nodes > 1000000U)
        {
            std::cout << "\rGenerating connections: " << counter << " / " << num_of_nodes;
        }
    }
    // std::cout << "\rGenerating connections: " << num_of_nodes << " / " << num_of_nodes << std::endl;
    return connections;
}