#include "infinite_adrcm_model.h"
#include "infinite_network.h"
#include "point.h"
#include "typedefs.h"

InfiniteAdrcmModel::InfiniteAdrcmModel(const py::array_t<double> &parameters_in, const uint32_t seed)
    : Model{seed}, InfiniteModel{}, AdrcmModel{parameters_in}
{
}

InfiniteNetwork InfiniteAdrcmModel::generate_network() const
{
    const auto vertices{create_vertices()};
    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    const auto connections{generate_connections(vertices)};
    const auto simplices{create_simplices_from_connections(connections)};
    const InfiniteNetwork network{max_dimension(), vertex_ids, simplices, 0};
    return network;
}

PointList InfiniteAdrcmModel::create_vertices() const
{
    std::uniform_real_distribution<float> uniform_distribution(0., 1.);

    const auto u{uniform_distribution(random_number_generator_)}; // birth time of the typical node

    // generate nodes
    PointList vertices{Point{0, u, 0.}};
    PointId id{1U};

    // generate older nodes which (u, 0) connects to
    // w is a random variable that is later transformed to be the birth time
    // there is a x0.5 due to parameter alpha, but there is a x2 as well due to +-
    const auto w_intensity_older_nodes{alpha() * beta() / (1. - gamma()) * std::pow(u, gamma() - 1.)};
    std::exponential_distribution<float> w_interarrival_time_distribution_older_nodes(w_intensity_older_nodes);
    auto w_older{w_interarrival_time_distribution_older_nodes(random_number_generator_)};
    while (w_older < std::pow(u, 1.F - gamma()))
    {
        // v: birth time of the neighbor of (u, 0)
        const auto v{std::pow(w_older, 1.F / (1.F - gamma()))};
        const auto radius{alpha() * beta() * std::pow(v, -gamma()) * std::pow(u, gamma() - 1.F)};

        // generate point uniformly
        Position position{radius * uniform_distribution(random_number_generator_)};

        // create point and save it
        vertices.push_back(Point{id, v, position});
        // increment w to arrive to the next birth time (before transformation)
        w_older += w_interarrival_time_distribution_older_nodes(random_number_generator_);
        ++id;
    }

    // generate younger nodes which connect to (u, 0)
    // w is a random variable that is later transformed to be the birth time
    const auto w_intensity_younger_nodes{alpha() * beta() / gamma() * std::pow(u, -gamma())};
    std::exponential_distribution<float> w_interarrival_time_distribution_younger_nodes{w_intensity_younger_nodes};
    auto w_younger{std::pow(u, gamma()) + w_interarrival_time_distribution_younger_nodes(random_number_generator_)};
    while (w_younger < 1.F)
    {
        // v: birth time of the neighbor of (u, 0)
        const auto v{std::pow(w_younger, 1.F / gamma())};
        const auto radius{alpha() * beta() * std::pow(u, -gamma()) * std::pow(v, gamma() - 1.F)};

        // generate point uniformly
        Position position(radius * uniform_distribution(random_number_generator_));

        // create point and save it
        vertices.push_back(Point{id, v, position});
        // increment w to arrive to the next birth time (before transformation)
        w_younger += w_interarrival_time_distribution_younger_nodes(random_number_generator_);
        ++id;
    }
    return vertices;
}
