
#include <vector>
#include <pybind11/numpy.h>

#include "tools.h"
#include "adrcm_model.h"
#include "network.h"

namespace py = pybind11;

AdrcmModel::Point::Point(const vertex_id id, const float birth_time, const float position)
    : id_(id), birth_time_(birth_time), position_(position)
{
}

AdrcmModel::Parameters::Parameters(const py::array_t<double> &model_parameters_input)
{
    const auto model_parameters = numpy_to_vector_1d<double>(model_parameters_input);
    max_dimension = static_cast<uint32_t>(model_parameters[0]);
    num_of_nodes = static_cast<uint32_t>(model_parameters[1]);
    alpha = model_parameters[2];
    beta = model_parameters[3];
    gamma = model_parameters[4];
    torus_dimension = static_cast<uint32_t>(model_parameters[5]);
}

AdrcmModel::AdrcmModel(
    const AdrcmModel::Parameters &parameters,
    const uint32_t seed) : parameters{parameters},
                           random_number_generator{seed}
{
}

Network AdrcmModel::generate_finite_network() const
{
    const auto vertices{create_vertices()};
    const auto connections{generate_network_connections_default(vertices, true)};
    const auto network{create_finite_network(vertices, connections)};
    return network;
}

std::vector<AdrcmModel::Point> AdrcmModel::create_vertices() const
{
    // generate nodes, nodes are added to this container in the order of their birth
    std::vector<AdrcmModel::Point> vertices;
    vertices.reserve(parameters.num_of_nodes);

    const auto birth_times{create_birth_times()};
    const auto time{birth_times.back()};

    std::uniform_real_distribution<> position_distribution{-time / 2., +time / 2.};
    for (auto id{0}; id < parameters.num_of_nodes; ++id)
    {
        const auto position{position_distribution(random_number_generator)};
        vertices.push_back(AdrcmModel::Point{
            id,
            birth_times.at(id) / time,
            static_cast<float>(position)});
    }
    return vertices;
}

std::vector<float> AdrcmModel::create_birth_times() const
{
    std::vector<float> birth_times;
    birth_times.reserve(parameters.num_of_nodes);

    std::exponential_distribution<> interarrival_time_distribution(1.);

    auto time{0.};
    for (auto id{0U}; id < parameters.num_of_nodes; ++id)
    {
        time += interarrival_time_distribution(random_number_generator);
        birth_times.push_back(time);
    }

    return birth_times;
}

connections AdrcmModel::generate_network_connections_default(
    const std::vector<AdrcmModel::Point> &vertices,
    const bool is_finite) const
{
    const auto torus_size{vertices.back().birth_time()};

    connections connections{};
    for (const auto &source : vertices)
    {
        for (auto target_node_id{0U}; target_node_id < source.id(); ++target_node_id)
        {
            const auto &target{vertices[target_node_id]};
            const auto birth_time_ratio{source.birth_time() / target.birth_time()};
            const auto distance{
                is_finite ? source.torus_distance(target, torus_size)
                          : source.distance(target)};
            const auto profile_function_argument{
                distance * source.birth_time() /
                (parameters.beta * pow(birth_time_ratio, parameters.gamma))};
            if (profile_function_argument < 0.5)
            {
                connections.push_back(std::pair(source.id(), target.id()));
            }
        }
    }
    return connections;
}

Network AdrcmModel::create_finite_network(const std::vector<AdrcmModel::Point> &vertices, const connections &connections) const
{
    std::vector<Simplex> simplices{};
    for (const auto &vertex : vertices)
    {
        simplices.push_back(Simplex(Vertex(vertex.id())));
    }
    for (const auto &connection : connections)
    {
        simplices.push_back(Simplex(Vertex(connection.first), Vertex(connection.second)));
    }

    const SimplicialComplex simplicial_complex{simplices.begin(), simplices.end(), true};
    const Network network{simplicial_complex};
    return network;
}

auto AdrcmModel::generate_infinite_networks(
    const uint32_t num_of_infinite_networks) const
{
    std::vector<connections> result{};
    const auto b{parameters.beta};
    const auto g{parameters.gamma};

    constexpr auto minimum_u_to_simulate{1e-7};
    std::uniform_real_distribution<> uniform_distribution_u(minimum_u_to_simulate, 1.);
    std::uniform_real_distribution<> uniform_distribution_y(-1., 1.);

    for (auto network_index{0U}; network_index < num_of_infinite_networks; ++network_index)
    {
        const auto u{uniform_distribution_u(random_number_generator)}; // birth time of the typical node

        // generate nodes, nodes are added to this container in the order of their birth
        std::vector<Point> nodes;

        // generate older nodes which (u, 0) connects to
        // w is a random variable that is later transformed to be the birth time
        // there is a x0.5 due to parameter alpha, but there is a x2 as well due to +-
        const auto w_intensity_older_nodes{b / (1. - g) * std::pow(u, g - 1.)};
        std::exponential_distribution<> w_interarrival_time_distribution_older_nodes(w_intensity_older_nodes);
        double w_older{w_interarrival_time_distribution_older_nodes(random_number_generator)};
        while (w_older < std::pow(u, 1. - g))
        {
            // v: birth time of the neighbor of (u, 0)
            const auto v{static_cast<float>(std::pow(w_older, 1. / (1. - g)))};
            // y: position of the neighbor of (u, 0)
            const auto y{static_cast<float>(
                0.5 * uniform_distribution_y(random_number_generator) *
                b * std::pow(v, -g) * std::pow(u, g - 1.))};
            // create point and save it
            nodes.push_back(Point{static_cast<vertex_id>(nodes.size()), v, y});
            // increment w to arrive to the next birth time (before transformation)
            w_older += w_interarrival_time_distribution_older_nodes(random_number_generator);
        }

        nodes.push_back(Point(nodes.size(), u, 0.));

        // generate younger nodes which connect to (u, 0)
        // w is a random variable that is later transformed to be the birth time
        const auto w_intensity_younger_nodes{b / g * std::pow(u, -g)};
        std::exponential_distribution<> w_interarrival_time_distribution_younger_nodes{w_intensity_younger_nodes};
        double w_younger{std::pow(u, g) + w_interarrival_time_distribution_younger_nodes(random_number_generator)};
        while (w_younger < 1.)
        {
            // v: birth time of the neighbor of (u, 0)
            const auto v{static_cast<float>(std::pow(w_younger, 1. / g))};
            // y: position of the neighbor of (u, 0)
            const auto y{static_cast<float>(
                0.5 * uniform_distribution_y(random_number_generator) *
                b * std::pow(u, -g) * std::pow(v, g - 1.))};
            // create point and save it
            nodes.push_back(Point{static_cast<vertex_id>(nodes.size()), v, y});
            // increment w to arrive to the next birth time (before transformation)
            w_younger += w_interarrival_time_distribution_younger_nodes(random_number_generator);
        }

        const auto connections{generate_network_connections_default(nodes, false)};
        result.push_back(connections);
    }

    return result;
}
