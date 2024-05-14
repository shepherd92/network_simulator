#include <cassert>
#include <mutex>

#include "infinite_network.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ISimplexList &nonempty_interactions,
    const uint32_t num_of_empty_interactions,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices, SimplexList{nonempty_interactions}, num_of_empty_interactions},
      typical_vertex_mark_{typical_vertex_mark_},
      marks_{marks},
      neighbors_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
}

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &nonempty_interactions,
    const uint32_t num_of_empty_interactions,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices, nonempty_interactions, num_of_empty_interactions},
      typical_vertex_mark_{typical_vertex_mark_},
      marks_{marks},
      neighbors_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
}

Mark InfiniteNetwork::typical_mark() const
{
    return typical_vertex_mark_;
}

SimplexList InfiniteNetwork::calc_simplices(const Dimension dimension)
{
    if (dimension == 0U)
    {
        return {};
    }
    const auto cofaces_of_typical_vertex{get_neighbors(dimension)};

    std::vector<Simplex> result{};
    for (auto &simplex : cofaces_of_typical_vertex)
    {
        const auto typical_vertex_is_oldest{
            std::all_of(
                std::execution::seq,
                simplex.vertices().begin(),
                simplex.vertices().end(),
                [&](const auto &vertex)
                {
                    return marks_[vertex] > typical_vertex_mark_;
                })};
        if (typical_vertex_is_oldest)
        {
            result.push_back(simplex);
        }
    }
    return result;
}

std::vector<uint32_t> InfiniteNetwork::calc_vertex_interaction_degree_distribution() const
{
    // assumption: all interactions are connected to the typical vertex
    return {nonempty_interactions_.size() + num_of_empty_interactions_};
}

const SimplexList &InfiniteNetwork::get_neighbors(const Dimension dimension)
{
    assert(dimension <= max_dimension_);
    if (!neighbors_[dimension].has_value())
    {
        neighbors_[dimension] = calc_neighbors(dimension);
    }
    return *neighbors_[dimension];
}

SimplexList InfiniteNetwork::calc_neighbors(const Dimension dimension)
{
    if (dimension == 0U)
    {
        return SimplexList{};
    }
    // -1 because the typical vertex is implicitly included in the neighbors
    const auto neighbors{get_facets().faces(dimension - 1U)};
    return neighbors;
}

std::vector<uint32_t> InfiniteNetwork::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(neighbor_dimension > simplex_dimension);
    if (simplex_dimension == 0U && neighbor_dimension == 1U)
    {
        return {num_vertices()};
    }

    const auto &possible_cofaces{get_neighbors(neighbor_dimension)};
    if (simplex_dimension == 0U)
    {
        return {possible_cofaces.size()};
    }
    return Network::calc_degree_sequence(simplex_dimension, neighbor_dimension);
}

InfiniteNetwork InfiniteNetwork::filter(const PointIdList &vertices) const
{
    InfiniteNetwork result{*this};
    result.keep_only_vertices(vertices);
    return result;
}
