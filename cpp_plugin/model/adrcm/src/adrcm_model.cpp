#include <atomic>

#include "adrcm_model.h"
#include "point.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

AdrcmModel::Parameters::Parameters(const std::vector<double> &model_parameters)
{
    max_dimension = static_cast<Dimension>(model_parameters[0]);
    network_size = model_parameters[1];
    alpha = model_parameters[2];
    beta = model_parameters[3];
    gamma = model_parameters[4];
}

AdrcmModel::AdrcmModel(const std::vector<double> &parameters)
    : parameters_{parameters}
{
}

ConnectionList AdrcmModel::generate_connections(const PointList &vertices) const
{
    // create aliases
    const auto num_of_nodes{vertices.size()};

    std::uniform_real_distribution<> uniform_distribution(0., 1.);
    ConnectionList connections{};
    std::atomic<uint32_t> counter{0U};
    for (auto target_id{0U}; target_id < num_of_nodes; ++target_id)
    {
        const auto &target_vertex{vertices[target_id]};
        const auto s{target_vertex.mark()};

        const auto size_of_neighborhood{(alpha() * beta() / s)};

        for (auto source_id{target_id + 1U}; source_id < num_of_nodes; ++source_id)
        {
            const auto &source_node{vertices[source_id]};
            const auto distance_{distance(source_node, target_vertex)};
            if (distance_ > size_of_neighborhood)
            {
                continue;
            }

            const auto t{source_node.mark()};
            const auto max_distance_of_connection{(alpha() * beta() / t) * std::pow(t / s, gamma())};

            if (distance_ < max_distance_of_connection)
            {
                if (is_default_alpha())
                {
                    connections.push_back(std::pair(source_id, target_id));
                }
                else
                {
                    const auto random_number{uniform_distribution(random_number_generator_)};
                    if (random_number < 1 / (2. * alpha()))
                    {
                        connections.push_back(std::pair(source_id, target_id));
                    }
                }
            }
        }
        log_progress(++counter, num_of_nodes, 10000U, "Generating connections");
    }
    log_progress(counter, num_of_nodes, 1U, "Generating connections");
    return connections;
}

SimplexList AdrcmModel::create_simplices_from_connections(const ConnectionList &connections) const
{
    std::vector<Simplex> simplices{};
    simplices.reserve(connections.size());
    for (const auto &connection : connections)
    {
        simplices.emplace_back(Simplex{PointIdList{connection.first, connection.second}});
    }
    return SimplexList{simplices};
}

Dimension AdrcmModel::max_dimension() const
{
    return parameters_.max_dimension;
}

float AdrcmModel::lambda() const
{
    return parameters_.network_size;
}

float AdrcmModel::alpha() const
{
    return parameters_.alpha;
}

float AdrcmModel::beta() const
{
    return parameters_.beta;
}

float AdrcmModel::gamma() const
{
    return parameters_.gamma;
}

bool AdrcmModel::is_default_alpha() const
{
    return is_close(alpha(), 0.5);
}