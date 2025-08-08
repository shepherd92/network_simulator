#include <cassert>
#include <mutex>

#include "infinite_hypergraph.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"

InfiniteHypergraph::InfiniteHypergraph(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &interactions,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices},
      InfiniteNetwork{typical_vertex_mark_, marks},
      Hypergraph{interactions}
{
}

std::vector<uint32_t> InfiniteHypergraph::calc_simplex_interaction_degree_sequence(
    const Dimension simplex_dimension)
{
    if (simplex_dimension == 0)
    {
        // assumption: all interactions are connected to the typical vertex
        return {interactions_.size()};
    }

    // -1 because the typical vertex is implicitly included
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

std::vector<uint32_t> InfiniteHypergraph::calc_vertex_interaction_degree_distribution() const
{
    // assumption: all interactions are connected to the typical vertex
    return {interactions_.size()};
}

SimplexList InfiniteHypergraph::calc_neighbors(const Dimension dimension)
{
    // typical vertex is implicitly included in the neighbors
    return dimension == 0U ? SimplexList{} : interactions_.faces(dimension - 1);
}

std::vector<Dimension> InfiniteHypergraph::calc_interaction_dimension_distribution() const
{
    // calculate the distribution where the typical vertex is not included
    auto result{Hypergraph::calc_interaction_dimension_distribution()};

    // increment each element by one to account for the typical vertex
    for (auto &dimension : result)
    {
        ++dimension;
    }

    return result;
}
