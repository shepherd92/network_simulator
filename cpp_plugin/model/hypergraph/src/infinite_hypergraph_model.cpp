#include <atomic>
#include <iostream>
#include <mutex>
#include <tuple>

#include "infinite_hypergraph_model.h"
#include "infinite_network.h"
#include "neighborhood.h"
#include "point.h"
#include "rectangle.h"
#include "tools.h"
#include "typedefs.h"

InfiniteHypergraphModel::InfiniteHypergraphModel(const std::vector<double> &parameters_in, const uint32_t seed)
    : Model{seed}, InfiniteModel{seed}, HypergraphModel{parameters_in}
{
}

std::vector<std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>>
InfiniteHypergraphModel::generate_networks(const uint32_t num_of_infinite_networks) const
{
    std::vector<std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>> result{};
    result.reserve(num_of_infinite_networks);

    for (auto network_index{0U}; network_index < num_of_infinite_networks; ++network_index)
    {
        const auto network_interface{generate_network()};
        result.push_back(network_interface);
        log_progress(network_index, num_of_infinite_networks, 1U, "Generating infinite networks");
    }
    return result;
}

std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>
InfiniteHypergraphModel::generate_network() const
{
    const PointId typical_vertex_id{0};
    std::uniform_real_distribution<Mark> mark_distribution(0., 1.);
    const auto u{mark_distribution(random_number_generator_)}; // mark of the typical node

    const auto interactions{create_interactions(u)};
    const auto interaction_mark_position_pairs{convert_to_mark_position_pairs(interactions)};
    const auto transformed_interactions{transform_interactions(interactions)};
    const std::vector<std::unique_ptr<NeighborhoodPart>> dominating_neighborhood_parts{
        create_dominating_neighborhood_parts(transformed_interactions)};

    auto vertices{create_vertices(interactions)};
    vertices.emplace_back(Point{u, 0.F, typical_vertex_id}); // add the typical node
    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};
    const auto vertex_marks{convert_to_mark_list(vertices)};

    const auto connections{generate_connections(vertices, interactions)};
    const auto simplices{create_simplices_from_connections(connections)};

    const InfiniteNetwork network{max_dimension(), vertex_ids, simplices, typical_vertex_id, vertex_marks};

    return {network, vertex_mark_position_pairs, interaction_mark_position_pairs};
}

PointList InfiniteHypergraphModel::create_interactions(const Mark u) const
{
    const PointId typical_vertex_id{0};
    // expected number of connecting interactions: 2 * b * l' * u^(-g) / (1 - g')
    std::poisson_distribution<int32_t> poisson_distribution_interactions(
        2 * beta() * lambda_prime() * std::pow(u, -gamma()) / (1 - gamma_prime()));
    const auto num_of_interactions{poisson_distribution_interactions(random_number_generator_)};
    const auto interaction_marks{generate_marks(num_of_interactions, MIN_MARK)};
    const auto interaction_positions{generate_positions_in_vertex_neighborhood(Point{u, 0.F, typical_vertex_id}, interaction_marks)};

    PointList interactions{};
    interactions.reserve(num_of_interactions);
    for (auto index{0}; index < num_of_interactions; ++index)
    {
        interactions.emplace_back(Point{interaction_marks[index], interaction_positions[index], index});
    }

    return interactions;
}

PointList InfiniteHypergraphModel::transform_interactions(const PointList &interactions) const
{
    PointList transformed_interactions{};
    transformed_interactions.reserve(interactions.size());
    std::transform(
        execution_policy,
        interactions.begin(), interactions.end(),
        std::back_inserter(transformed_interactions),
        [&](const auto &point)
        {
            return Point{beta() * std::pow(point.mark(), -gamma_prime()), point.position(), point.id()};
        });
    return transformed_interactions;
}

std::vector<std::unique_ptr<NeighborhoodPart>>
InfiniteHypergraphModel::create_dominating_neighborhood_parts(const PointList &transformed_interactions) const
{
    auto centers{create_neighborhood_centers(transformed_interactions)};
    merge_neighborhood_centers(centers);
    auto tails{create_neighborhood_tails(transformed_interactions)};
    remove_center_domains_from_neighborhood_tails(tails, centers);
    const auto dominating_tails{determine_dominating_hyperbolas(tails)};

    std::vector<std::unique_ptr<NeighborhoodPart>> dominating_neighborhood_parts{};
    dominating_neighborhood_parts.reserve(centers.size() + dominating_tails.size());
    std::transform(
        execution_policy,
        centers.begin(), centers.end(),
        std::back_inserter(dominating_neighborhood_parts),
        [](const auto &center)
        {
            return std::make_unique<NeighborhoodPart>(center);
        });
    std::transform(
        execution_policy,
        dominating_tails.begin(), dominating_tails.end(),
        std::back_inserter(dominating_neighborhood_parts),
        [](const auto &tail)
        {
            return std::make_unique<NeighborhoodPart>(tail);
        });
    return dominating_neighborhood_parts;
}

std::vector<Center>
InfiniteHypergraphModel::create_neighborhood_centers(const PointList &transformed_interactions) const
{
    std::vector<Center> centers{};
    centers.reserve(transformed_interactions.size());
    std::transform(
        execution_policy,
        transformed_interactions.begin(), transformed_interactions.end(),
        std::back_inserter(centers),
        [&](const auto &interaction)
        {
            return get_neighborhood_center(interaction);
        });
}

void InfiniteHypergraphModel::merge_neighborhood_centers(std::vector<Center> &centers) const
{
    std::sort( // sort by left border
        execution_policy,
        centers.begin(), centers.end(),
        [](const auto &center1, const auto &center2)
        {
            return center1.left() < center2.left();
        });

    std::vector<Center> merged_centers{centers[0]};
    for (auto index{1U}; index < centers.size(); ++index)
    {
        if (centers[index].left() <= merged_centers.back().right())
        {
            merged_centers.back().set_right(centers[index].right());
        }
        else
        {
            merged_centers.push_back(centers[index]);
        }
    }
    centers = merged_centers;
}

std::vector<Hyperbola>
InfiniteHypergraphModel::create_neighborhood_tails(const PointList &transformed_interactions) const
{
    std::vector<Hyperbola> tails{};
    tails.reserve(2 * transformed_interactions.size());

    std::for_each(
        std::execution::seq,
        transformed_interactions.begin(), transformed_interactions.end(),
        [&](const auto &interaction)
        {
            tails.push_back(get_neighborhood_left_tail(interaction));
            tails.push_back(get_neighborhood_right_tail(interaction));
        });
}

void InfiniteHypergraphModel::remove_center_domains_from_neighborhood_tails(std::vector<Hyperbola> &tails, const std::vector<Center> &centers) const
{
    std::vector<Hyperbola> result{};
    std::for_each(
        std::execution::seq,
        tails.begin(), tails.end(),
        [&](auto &tail)
        {
            std::vector<Hyperbola> remaining_parts{};
            for (const auto &center : centers)
            {
                if (tail.right() < center.left())
                {
                    // center:         |-----|
                    // tail:   |-----|
                    remaining_parts.push_back(tail);
                    // we reached the right end of the tail
                    break;
                }
                else if (tail.left() < center.left() && tail.right() < center.right())
                {
                    // center:    |-----|
                    // tail:   |-----|
                    remaining_parts.push_back(Hyperbola(tail.left(), center.left(), tail.position(), tail.transformed_mark()));
                    // we reached the right end of the tail
                    break;
                }
                else if (tail.left() < center.left() && tail.right() > center.right())
                {
                    // center:    |-----|
                    // tail:   |-----------|
                    remaining_parts.push_back(Hyperbola(tail.left(), center.left(), tail.position(), tail.transformed_mark()));
                    tail.set_left(center.right());
                }
                else if (center.left() < tail.left() && center.right() > tail.right())
                {
                    // center: |-------|
                    // tail:     |---|
                    // we reached the right end of the tail
                    break;
                }
                else if (center.right() < tail.left())
                {
                    // center: |-----|
                    // tail:            |-----|
                    continue; // do nothing
                }
                else if (center.left() < tail.left() && center.right() < tail.right())
                {
                    // center: |-----|
                    // tail:      |-----|
                    tail.set_left(center.right());
                }
                else
                {
                    assert(false);
                }
            }
            if (tail.side() == Hyperbola::Side::RIGHT)
            {
                // if the tail is right sided, then put the last part to the remaining parts
                remaining_parts.push_back(tail);
            }
            result.insert(result.end(), remaining_parts.begin(), remaining_parts.end());
        });
    tails = result;
}

std::vector<Hyperbola>
InfiniteHypergraphModel::determine_dominating_hyperbolas(std::vector<Hyperbola> &tails) const
{
    std::sort(
        execution_policy,
        tails.begin(), tails.end(),
        [](const auto &tail1, const auto &tail2)
        {
            return tail1.left() < tail2.left();
        });

    std::vector<Hyperbola> result{};
    Position currently_left{tails[0].left()};
    for (auto index{0U}; index < tails.size() - 1; ++index)
    {
        std::vector<Hyperbola> hyperbolas_in_this_domain{};
        if (is_close(tails[index].left(), currently_left))
        {
            hyperbolas_in_this_domain.push_back(tails[index]);
        }
        else
        {
            // all hyperbolas in the domain are collected
            // PROCESSING
            currently_left = tails[index].left();
            hyperbolas_in_this_domain.clear();
            hyperbolas_in_this_domain.push_back(tails[index]);
        }
    }
    return result;
}

// PointList InfiniteHypergraphModel::create_vertices(const PointList &interactions) const
// {
//     PointList vertices{};
//     std::mutex mutex{};
//     std::atomic<PointId> vertex_id{1U}; // 0 is reserved for the typical vertex
//     for (auto interaction_index{0U}; interaction_index < interactions.size(); ++interaction_index)
//     {
//         const auto vertices_in_neighborhood{create_vertices_in_interaction_neighborhood(interactions[interaction_index])};
//         std::for_each(
//             execution_policy,
//             vertices_in_neighborhood.begin(), vertices_in_neighborhood.end(),
//             [&](const auto &potential_vertex)
//             {
//                 auto should_be_discarded{false};
//                 for (auto index{0U}; index < interaction_index; ++index)
//                 {
//                     if (connects(potential_vertex, interactions[index]))
//                     {
//                         std::lock_guard<std::mutex> lock{mutex};
//                         should_be_discarded = true;
//                         break;
//                     }
//                 }
//                 if (!should_be_discarded)
//                 {
//                     std::lock_guard<std::mutex> lock{mutex};
//                     vertices.emplace_back(Point(potential_vertex.mark(), potential_vertex.position(), vertex_id));
//                     ++vertex_id;
//                 }
//             });
//     }
//     return vertices;
// }

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

Hyperbola InfiniteHypergraphModel::get_neighborhood_left_tail(const Point &interaction) const
{
    return Hyperbola(
        -std::numeric_limits<float>::infinity(),
        interaction.position() - interaction.mark(),
        interaction.position(),
        interaction.mark());
}

Center InfiniteHypergraphModel::get_neighborhood_center(const Point &interaction) const
{
    return Center(
        interaction.position() - interaction.mark(),
        interaction.position() + interaction.mark());
}

Hyperbola InfiniteHypergraphModel::get_neighborhood_right_tail(const Point &interaction) const
{
    return Hyperbola(
        interaction.position() + interaction.mark(),
        std::numeric_limits<float>::infinity(),
        interaction.position(),
        interaction.mark());
}
