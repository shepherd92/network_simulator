#include <cassert>
#include <mutex>

#include "infinite_network.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"

InfiniteNetwork::InfiniteNetwork(
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : typical_vertex_mark_{typical_vertex_mark_},
      marks_{marks},
      neighbors_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
}

InfiniteNetwork::InfiniteNetwork(InfiniteNetwork &&other) noexcept
{
    typical_vertex_mark_ = std::move(other.typical_vertex_mark_);
    marks_ = std::move(other.marks_);
    neighbors_ = std::move(other.neighbors_);
}

InfiniteNetwork &InfiniteNetwork::operator=(InfiniteNetwork &&other) noexcept
{
    if (this != &other)
    {
        typical_vertex_mark_ = std::move(other.typical_vertex_mark_);
        marks_ = std::move(other.marks_);
        neighbors_ = std::move(other.neighbors_);
    }
    return *this;
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
    // This method calculates those simplices where the typical vertex is the lowest-mark vertex.
    // The result contains all these simplices without the typical vertex.
    if (dimension == 0U)
    {
        return {}; // typical vertex is not included in the simplices
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

std::vector<uint32_t> InfiniteNetwork::calc_coface_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(neighbor_dimension > simplex_dimension);
    if (simplex_dimension == 0U && neighbor_dimension == 1U)
    {
        // return number of all vertices, not just the ones younger than the typical vertex
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
