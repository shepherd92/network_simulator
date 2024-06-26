#include <atomic>
#include <cassert>
#include <iostream>
#include <mutex>
#include <thread>
#include <tuple>

#include "infinite_hypergraph_model.h"
#include "infinite_network.h"
#include "neighborhood.h"
#include "point.h"
#include "rectangle.h"
#include "tools.h"
#include "typedefs.h"

InfiniteHypergraphModel::InfiniteHypergraphModel(const std::vector<double> &parameters_in, const uint32_t seed)
    : Model{seed}, InfiniteModel{seed}, HypergraphModel{parameters_in}, interactions_only_{parameters_in[7] > 0.5}
{
}

std::vector<std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>>
InfiniteHypergraphModel::generate_networks(const uint32_t num_of_infinite_networks) const
{
    std::uniform_real_distribution<Mark> mark_distribution(0.F, MAX_MARK);

    std::vector<std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>> result{};
    result.reserve(num_of_infinite_networks);

    for (auto network_index{0U}; network_index < num_of_infinite_networks; ++network_index)
    {
        const auto typical_vertex_mark{std::max(mark_distribution(random_number_generator_), MIN_MARK)}; // mark of the typical node
        const auto network_interface{generate_network(typical_vertex_mark)};
        result.push_back(network_interface);
        log_progress(network_index, num_of_infinite_networks, 1U, "Generating infinite networks");
    }
    return result;
}

std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>
InfiniteHypergraphModel::generate_network(const Mark typical_vertex_mark) const
{
    const auto interactions{create_interactions(typical_vertex_mark)};
    const auto interaction_ids{convert_to_id_list(interactions)};
    const auto interaction_mark_position_pairs{convert_to_mark_position_pairs(interactions)};

    if (interactions_only())
    {
        std::vector<Simplex> simplices{interactions.size(), PointIdList{}};
        return {InfiniteNetwork{max_dimension(), {}, SimplexList{simplices}, typical_vertex_mark, {}}, {}, interaction_mark_position_pairs};
    }

    PointList vertices{};
    if (!interactions.empty())
    {
        const auto transformed_interactions{transform_interactions(interactions)};
        const auto dominating_neighborhood_parts{create_dominating_neighborhood_parts(transformed_interactions)};
        const auto ordinary_vertices{create_vertices(dominating_neighborhood_parts.first, dominating_neighborhood_parts.second)};
        vertices = std::move(ordinary_vertices);
    }

    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto vertex_marks{convert_to_mark_list(vertices)};
    const auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    const auto connections{generate_connections(vertices, interactions)};
    const auto simplices{create_interaction_simplices_from_connections(interaction_ids, connections)};

    const InfiniteNetwork network{max_dimension(), vertex_ids, simplices, typical_vertex_mark, vertex_marks};

    return {network, vertex_mark_position_pairs, interaction_mark_position_pairs};

    // check if the connections are correct
    // for (const auto &vertex : vertices)
    // {
    //     for (const auto &interaction : interactions)
    //     {
    //         const bool should_connect{
    //             std::abs(vertex.position() - interaction.position()) < beta() * std::pow(vertex.mark(), -gamma()) * std::pow(interaction.mark(), -gamma_prime())};
    //         const bool in_connections{
    //             std::find_if(
    //                 connections.begin(), connections.end(),
    //                 [&](const auto &connection)
    //                 { return connection.first == vertex.id() && connection.second == interaction.id(); }) != connections.end()};
    //         assert(should_connect == in_connections);
    //     }
    // }
}

PointList InfiniteHypergraphModel::create_interactions(const Mark typical_vertex_mark) const
{
    // expected number of connecting interactions: 2 * b * l' * u^(-g) / (1 - g')
    std::poisson_distribution<int32_t> poisson_distribution_interactions(
        2 * beta() * lambda_prime() * std::pow(typical_vertex_mark, -gamma()) / (1 - gamma_prime()));
    const auto num_of_interactions{poisson_distribution_interactions(random_number_generator_)};
    const auto interaction_marks{generate_marks(num_of_interactions, MIN_MARK)};
    const auto interaction_positions{generate_positions_in_vertex_neighborhood(Point{typical_vertex_mark, 0.F}, interaction_marks)};

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
    std::for_each(
        std::execution::seq,
        interactions.begin(), interactions.end(),
        [&](const auto &point)
        {
            transformed_interactions.emplace_back(Point{beta() * std::pow(point.mark(), -gamma_prime()), point.position()});
        });
    return transformed_interactions;
}

std::pair<std::vector<Center>, std::vector<Hyperbola>>
InfiniteHypergraphModel::create_dominating_neighborhood_parts(const PointList &transformed_interactions) const
{
    if (transformed_interactions.empty())
    {
        return {};
    }
    auto centers{create_neighborhood_centers(transformed_interactions)};
    auto tails{create_neighborhood_tails(transformed_interactions)};
    const auto tail_parts_in_slots{remove_center_domains_from_neighborhood_tails(tails, centers)};
    const auto dominating_tails{determine_dominating_hyperbolas(tail_parts_in_slots)};
    return std::make_pair(centers, dominating_tails);
}

std::vector<Center>
InfiniteHypergraphModel::create_neighborhood_centers(const PointList &transformed_interactions) const
{
    std::vector<Center> centers{};
    centers.reserve(transformed_interactions.size());
    std::for_each(
        std::execution::seq,
        transformed_interactions.begin(), transformed_interactions.end(),
        [&](const auto &interaction)
        {
            centers.emplace_back(get_neighborhood_center(interaction));
        });
    return merge_neighborhood_centers(centers);
}

std::vector<Center> InfiniteHypergraphModel::merge_neighborhood_centers(std::vector<Center> &centers) const
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
        if (is_significantly_less(centers[index].left(), merged_centers.back().right()))
        {
            merged_centers.back().set_right(centers[index].right());
        }
        else
        {
            merged_centers.push_back(centers[index]);
        }
    }
    return merged_centers;
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
    return tails;
}

std::vector<std::vector<Hyperbola>>
InfiniteHypergraphModel::remove_center_domains_from_neighborhood_tails(
    const std::vector<Hyperbola> &tails,
    const std::vector<Center> &centers) const
{
    /*Algorithm:
    Assumption: centers are disjoint and sorted
    Extract hyperbolas in each slot (centers: |---|)
       slot1 |-------|  slot2  |--|  slot3  |---------| slot4 |----|  slot5*/
    std::vector<std::vector<Hyperbola>> result{};
    Position current_slot_left{-std::numeric_limits<float>::infinity()};
    std::for_each(
        std::execution::seq,
        centers.begin(), centers.end(),
        [&](const auto &center)
        {
            std::vector<Hyperbola> hyperbolas_this_slot{};
            Position current_slot_right{center.left()};
            std::for_each(
                std::execution::seq,
                tails.begin(), tails.end(),
                [&](auto &tail)
                {
                    const auto tail_part{tail.intersect_domain(current_slot_left, current_slot_right)};
                    if (tail_part.has_value())
                    {
                        hyperbolas_this_slot.push_back(tail_part.value());
                    }
                });
            current_slot_left = center.right();
            result.push_back(hyperbolas_this_slot);
        });

    // add the last slot from <last center>.right() to infinity
    std::vector<Hyperbola> hyperbolas_last_slot{};
    std::for_each(
        std::execution::seq,
        tails.begin(), tails.end(),
        [&](auto &tail)
        {
            const auto tail_part{tail.intersect_domain(current_slot_left, std::numeric_limits<float>::infinity())};
            if (tail_part.has_value())
            {
                hyperbolas_last_slot.push_back(tail_part.value());
            }
        });
    result.push_back(hyperbolas_last_slot);
    return result;
}

std::vector<Hyperbola>
InfiniteHypergraphModel::determine_dominating_hyperbolas(
    const std::vector<std::vector<Hyperbola>> &tails_in_slots) const
{
    std::vector<Hyperbola> result{};
    for (const auto &slot : tails_in_slots)
    {
        std::vector<Hyperbola> dominating_hyperbolas{slot[0]};

        for (const auto &hyperbola : slot)
        {
            std::mutex mutex{};
            std::vector<Hyperbola> new_dominating_hyperbolas{};
            std::for_each(
                execution_policy,
                dominating_hyperbolas.begin(), dominating_hyperbolas.end(),
                [&](const auto &dominating_hyperbola)
                {
                    const auto dominating_hyperbola_parts{dominating_hyperbola.get_dominating_hyperbola_parts(hyperbola)};
                    std::lock_guard<std::mutex> lock(mutex);
                    new_dominating_hyperbolas.insert(
                        new_dominating_hyperbolas.end(),
                        dominating_hyperbola_parts.begin(),
                        dominating_hyperbola_parts.end());
                });
            dominating_hyperbolas = new_dominating_hyperbolas;
        }
        result.insert(result.end(), dominating_hyperbolas.begin(), dominating_hyperbolas.end());
    }
    return result;
}

PointList InfiniteHypergraphModel::create_vertices(
    const std::vector<Center> &centers,
    const std::vector<Hyperbola> &hyperbolas) const
{
    PointList vertices{};
    PointList vertices_under_centers{create_vertices_under_centers(centers)};
    vertices.insert(vertices.end(), vertices_under_centers.begin(), vertices_under_centers.end());
    PointList vertices_under_tails{create_vertices_under_tails(hyperbolas)};
    vertices.insert(vertices.end(), vertices_under_tails.begin(), vertices_under_tails.end());

    PointId vertex_id{0U};
    for (auto &vertex : vertices)
    {
        vertex.set_id(vertex_id);
        ++vertex_id;
    }
    return vertices;
}

PointList InfiniteHypergraphModel::create_vertices_under_centers(const std::vector<Center> &centers) const
{
    PointList vertices{};
    for (const auto &center : centers)
    {
        auto new_vertices{center.create_points(parameters(), random_number_generator_)};
        vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
    }
    return vertices;
}

PointList InfiniteHypergraphModel::create_vertices_under_tails(const std::vector<Hyperbola> &hyperbolas) const
{
    PointList vertices{};
    for (const auto &hyperbola : hyperbolas)
    {
        auto new_vertices{hyperbola.create_points(parameters(), random_number_generator_)};
        vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
    }
    return vertices;
}

PositionList InfiniteHypergraphModel::generate_positions_in_vertex_neighborhood(
    const Point &vertex,
    const MarkList &marks) const
{
    std::uniform_real_distribution<Position> uniform_distribution(-1., 1.);
    const auto beta_x_mark_to_gamma{beta() * std::pow(vertex.mark(), -gamma())};

    PositionList positions{};
    positions.reserve(marks.size());
    std::for_each(
        std::execution::seq,
        marks.begin(), marks.end(),
        [&](const auto mark)
        {
            const auto relative_position{uniform_distribution(random_number_generator_) * beta_x_mark_to_gamma * std::pow(mark, -gamma_prime())};
            positions.push_back(vertex.position() + relative_position);
        });
    return positions;
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

bool InfiniteHypergraphModel::interactions_only() const
{
    return interactions_only_;
}

std::vector<uint32_t> InfiniteHypergraphModel::calc_vertex_interaction_degree_sequence_directly(const uint32_t num_of_networks)
{
    const float constant_term{2.F * beta() * lambda_prime() / (1.F - gamma_prime())};
    // std::uniform_real_distribution<Mark> mark_distribution(0.F, 1.F);
    // std::poisson_distribution<uint32_t> poisson_distribution(0.F);

    std::vector<uint32_t> result{};
    result.reserve(num_of_networks);
    std::uniform_real_distribution<Mark> mark_distribution(0.F, 1.F);

    for (auto i{0U}; i < num_of_networks; ++i)
    {
        const Mark typical_mark{mark_distribution(random_number_generator_)};
        const auto poisson_lambda{constant_term * fast_power(typical_mark, -gamma())};
        std::poisson_distribution<uint32_t> poisson_distribution(poisson_lambda);
        result.push_back(poisson_distribution(random_number_generator_));
    }
    return result;
}
