#include <cassert>
#include <mutex>

#include "infinite_network.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ISimplexList &interactions,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices, SimplexList{interactions}},
      typical_vertex_mark_{typical_vertex_mark_},
      marks_{marks},
      neighbors_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
}

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &interactions,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices, interactions},
      typical_vertex_mark_{typical_vertex_mark_},
      marks_{marks},
      neighbors_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
}

Mark InfiniteNetwork::typical_mark() const
{
    return typical_vertex_mark_;
}

SimplexList InfiniteNetwork::get_skeleton(const Dimension max_dimension)
{
    SimplexList result{};
    for (auto dimension{0}; dimension <= max_dimension; ++dimension)
    {
        result += get_neighbors(dimension);
    }
    return result;
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

std::vector<uint32_t> InfiniteNetwork::calc_simplex_interaction_degree_sequence(
    const Dimension simplex_dimension)
{
    if (simplex_dimension == 0)
    {
        // assumption: all interactions are connected to the typical vertex
        return {interactions_.size()};
    }

    auto simplex_degree_map{interactions_.calc_degree_sequence(simplex_dimension - 1)};

    // order of the degree values does not matter
    std::vector<uint32_t> result{};
    const auto &simplices{get_simplices(simplex_dimension)};
    result.reserve(simplex_degree_map.size());
    for (const auto &simplex : simplices)
    {
        // simplex has typical vertex, simplex_degree_map is for simplices without typical vertex
        result.emplace_back(simplex_degree_map[simplex]);
    }

    return result;
}

std::vector<uint32_t> InfiniteNetwork::calc_vertex_interaction_degree_distribution() const
{
    // assumption: all interactions are connected to the typical vertex
    return {};
}

std::vector<Dimension> InfiniteNetwork::calc_facet_dimension_distribution()
{
    // first calculate the distribution where the typical vertex is not included
    auto result{Network::calc_facet_dimension_distribution()};
    // add dimension 0 by shifting the vector (there are no facets of dimension 0)
    result.insert(result.begin(), 0U);
    return result;
}

std::vector<Dimension> InfiniteNetwork::calc_interaction_dimension_distribution() const
{
    // calculate the distribution where the typical vertex is not included
    auto result{Network::calc_interaction_dimension_distribution()};
    return result;
}

std::vector<uint32_t> InfiniteNetwork::calc_coface_degree_sequence(
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

    // -1 because the typical vertex is implicitly included
    auto simplex_degree_map{possible_cofaces.calc_degree_sequence(simplex_dimension - 1U)};
    std::vector<uint32_t> result{};
    const auto &simplices{get_simplices(simplex_dimension)};
    for (const auto &simplex : simplices)
    {
        result.emplace_back(simplex_degree_map[simplex]);
    }
    return result;
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

InfiniteNetwork InfiniteNetwork::filter(const PointIdList &vertices) const
{
    InfiniteNetwork result{*this};
    result.keep_only_vertices(vertices);
    return result;
}
