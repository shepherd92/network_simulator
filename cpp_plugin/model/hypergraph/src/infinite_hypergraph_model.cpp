#include <atomic>
#include <iostream>
#include <mutex>
#include <tuple>

#include "infinite_hypergraph_model.h"
#include "infinite_network.h"
#include "point.h"
#include "rectangle.h"
#include "typedefs.h"

InfiniteHypergraphModel::InfiniteHypergraphModel(const std::vector<double> &parameters_in, const uint32_t seed)
    : Model{seed}, HypergraphModel{parameters_in}
{
}

InfiniteNetwork InfiniteHypergraphModel::generate_network() const
{
    const PointId typical_vertex_id{0};
    std::uniform_real_distribution<Mark> mark_distribution(0., 1.);
    const auto u{mark_distribution(random_number_generator_)}; // mark of the typical node

    const auto interactions{create_interactions(u)};
    const auto interaction_mark_position_pairs{convert_to_mark_position_pairs(interactions)};

    auto vertices{create_vertices(interactions)};
    vertices.emplace_back(Point{0.F, u, typical_vertex_id}); // add the typical node
    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    const auto connections{generate_connections(vertices, interactions)};
    const auto simplices{create_simplices_from_connections(connections)};

    const InfiniteNetwork network{max_dimension(), vertex_ids, simplices, typical_vertex_id};

    return network;
}

PointList InfiniteHypergraphModel::create_interactions(const Mark u) const
{
    const PointId typical_vertex_id{0};
    // expected number of connecting interactions: 2 * b * l' * u^(-g) / (1 - g')
    std::poisson_distribution<int32_t> poisson_distribution_interactions(
        2 * beta() * lambda_prime() * std::pow(u, -gamma()) / (1 - gamma_prime()));
    const auto num_of_interactions{poisson_distribution_interactions(random_number_generator_)};
    const auto interaction_marks{generate_marks(num_of_interactions, MIN_MARK)};
    const auto interaction_positions{generate_positions_in_vertex_neighborhood(Point{0.F, u, typical_vertex_id}, interaction_marks)};

    PointList interactions{};
    interactions.reserve(num_of_interactions);
    for (auto index{0}; index < num_of_interactions; ++index)
    {
        interactions.emplace_back(Point{interaction_marks[index], interaction_positions[index], index});
    }

    return interactions;
}

PointList InfiniteHypergraphModel::create_vertices(const PointList &interactions) const
{
    PointList vertices{};
    std::mutex mutex{};
    std::atomic<PointId> vertex_id{1U}; // 0 is reserved for the typical vertex
    for (auto interaction_index{0U}; interaction_index < interactions.size(); ++interaction_index)
    {
        const auto vertices_in_neighborhood{create_vertices_in_interaction_neighborhood(interactions[interaction_index])};
        std::for_each(
            execution_policy,
            vertices_in_neighborhood.begin(), vertices_in_neighborhood.end(),
            [&](const auto &potential_vertex)
            {
                auto should_be_discarded{false};
                for (auto index{0U}; index < interaction_index; ++index)
                {
                    if (connects(potential_vertex, interactions[index]))
                    {
                        std::lock_guard<std::mutex> lock{mutex};
                        should_be_discarded = true;
                        break;
                    }
                }
                if (!should_be_discarded)
                {
                    std::lock_guard<std::mutex> lock{mutex};
                    vertices.emplace_back(Point(potential_vertex.mark(), potential_vertex.position(), vertex_id));
                    ++vertex_id;
                }
            });
    }
    return vertices;
}

PointList InfiniteHypergraphModel::create_vertices_in_interaction_neighborhood(const Point &interaction) const
{
    const auto expected_num_of_vertices{2 * beta() * lambda() * std::pow(interaction.mark(), -gamma_prime()) / (1. - gamma())};
    const auto num_of_vertices{std::poisson_distribution<int32_t>(expected_num_of_vertices)(random_number_generator_)};
    const auto marks{generate_marks(num_of_vertices, MIN_MARK)};
    const auto positions{generate_positions_in_interaction_neighborhood(interaction, marks)};
    PointList vertices{};
    vertices.reserve(num_of_vertices);
    for (auto index{0}; index < num_of_vertices; ++index)
    {
        vertices.emplace_back(Point{marks[index], positions[index], index});
    }
    return vertices;
}

PositionList InfiniteHypergraphModel::generate_positions_in_vertex_neighborhood(
    const Point &vertex,
    const MarkList &marks) const
{
    return generate_positions_in_neighborhood(vertex, marks, gamma(), gamma_prime());
}

PositionList InfiniteHypergraphModel::generate_positions_in_interaction_neighborhood(
    const Point &interaction,
    const MarkList &marks) const
{
    return generate_positions_in_neighborhood(interaction, marks, gamma_prime(), gamma());
}

PositionList InfiniteHypergraphModel::generate_positions_in_neighborhood(
    const Point &point,
    const MarkList &marks,
    const float exponent_of_central_point,
    const float exponent_of_points_in_neighborhood) const
{
    std::uniform_real_distribution<Position> uniform_distribution(-1., 1.);
    const auto beta_x_mark_to_gamma{beta() * std::pow(point.mark(), -exponent_of_central_point)};

    PositionList positions{};
    positions.reserve(marks.size());
    std::for_each(
        std::execution::seq,
        marks.begin(), marks.end(),
        [&](const auto mark)
        {
            const auto position{uniform_distribution(random_number_generator_) * beta_x_mark_to_gamma * std::pow(mark, -exponent_of_points_in_neighborhood)};
            positions.push_back(position);
        });
    return positions;
}

bool InfiniteHypergraphModel::rectangle_points_surely_connect(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const
{
    const Point vtl{vertex_rectangle.top(), vertex_rectangle.left()};
    const Point vtr{vertex_rectangle.top(), vertex_rectangle.right()};
    const Point itl{interaction_rectangle.top(), interaction_rectangle.left()};
    const Point itr{interaction_rectangle.top(), interaction_rectangle.right()};

    if (fabs(vtl.position() - itl.position()) < vtl.mark() * itl.mark() &&
        fabs(vtl.position() - itr.position()) < vtl.mark() * itr.mark() &&
        fabs(vtr.position() - itl.position()) < vtr.mark() * itl.mark() &&
        fabs(vtr.position() - itr.position()) < vtr.mark() * itr.mark())
    {
        return true;
    }

    return false;
}