#include "infinite_adrcm_model.h"
#include "infinite_clique_complex.h"
#include "point.h"
#include "typedefs.h"

InfiniteAdrcmModel::InfiniteAdrcmModel(const std::vector<double> &parameters_in, const uint32_t seed)
    : Model{seed}, InfiniteModel{seed}, AdrcmModel{parameters_in}
{
}

std::vector<InfiniteCliqueComplex> InfiniteAdrcmModel::generate_networks(const uint32_t num_of_infinite_networks) const
{
    std::vector<InfiniteCliqueComplex> networks;
    networks.reserve(num_of_infinite_networks);
    for (uint32_t i = 0; i < num_of_infinite_networks; ++i)
    {
        networks.emplace_back(generate_network());
    }
    return networks;
}

InfiniteCliqueComplex InfiniteAdrcmModel::generate_network() const
{
    std::uniform_real_distribution<float> uniform_distribution(0., 1.);
    const auto typical_vertex_mark{uniform_distribution(random_number_generator_)}; // birth time of the typical node

    const auto vertices{create_vertices(typical_vertex_mark)};
    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};
    const auto vertex_marks{convert_to_mark_list(vertices)};

    const auto connections{generate_connections(vertices)};
    const InfiniteCliqueComplex network{max_dimension(), vertex_ids, typical_vertex_mark, vertex_marks};
    return network;
}

PointList InfiniteAdrcmModel::create_vertices(const Mark typical_vertex_mark) const
{
    std::uniform_real_distribution<float> uniform_distribution(0., 1.);
    // generate nodes
    PointList vertices{Point{typical_vertex_mark, 0.F, 0}};
    PointId id{1};

    // generate older nodes which (u, 0) connects to
    // w is a random variable that is later transformed to be the birth time
    // there is a x0.5 due to parameter alpha, but there is a x2 as well due to +-
    const auto w_intensity_older_nodes{alpha() * beta() / (1. - gamma()) * std::pow(typical_vertex_mark, gamma() - 1.)};
    std::exponential_distribution<float> w_interarrival_time_distribution_older_nodes(w_intensity_older_nodes);
    auto w_older{w_interarrival_time_distribution_older_nodes(random_number_generator_)};
    while (w_older < std::pow(typical_vertex_mark, 1.F - gamma()))
    {
        // v: birth time of the neighbor of (u, 0)
        const auto v{std::pow(w_older, 1.F / (1.F - gamma()))};
        const auto radius{alpha() * beta() * std::pow(v, -gamma()) * std::pow(typical_vertex_mark, gamma() - 1.F)};

        // generate point uniformly
        Position position{radius * uniform_distribution(random_number_generator_)};

        // create point and save it
        vertices.push_back(Point{v, position, id});
        // increment w to arrive to the next birth time (before transformation)
        w_older += w_interarrival_time_distribution_older_nodes(random_number_generator_);
        ++id;
    }

    // generate younger nodes which connect to (u, 0)
    // w is a random variable that is later transformed to be the birth time
    const auto w_intensity_younger_nodes{alpha() * beta() / gamma() * std::pow(typical_vertex_mark, -gamma())};
    std::exponential_distribution<float> w_interarrival_time_distribution_younger_nodes{w_intensity_younger_nodes};
    auto w_younger{std::pow(typical_vertex_mark, gamma()) + w_interarrival_time_distribution_younger_nodes(random_number_generator_)};
    while (w_younger < 1.F)
    {
        // v: birth time of the neighbor of (u, 0)
        const auto v{std::pow(w_younger, 1.F / gamma())};
        const auto radius{alpha() * beta() * std::pow(typical_vertex_mark, -gamma()) * std::pow(v, gamma() - 1.F)};

        // generate point uniformly
        Position position(radius * uniform_distribution(random_number_generator_));

        // create point and save it
        vertices.push_back(Point{v, position, id});
        // increment w to arrive to the next birth time (before transformation)
        w_younger += w_interarrival_time_distribution_younger_nodes(random_number_generator_);
        ++id;
    }
    return vertices;
}
