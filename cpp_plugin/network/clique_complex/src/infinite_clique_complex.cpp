#include <cassert>
#include <mutex>

#include "infinite_clique_complex.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"

InfiniteCliqueComplex::InfiniteCliqueComplex(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ConnectionList &edges,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices},
      InfiniteNetwork{typical_vertex_mark_, marks},
      CliqueComplex{edges}
{
}

InfiniteCliqueComplex::InfiniteCliqueComplex(InfiniteCliqueComplex &&other) noexcept
    : Network{std::move(other)},
      InfiniteNetwork{std::move(other)},
      CliqueComplex{std::move(other)}
{
}

// move assignment defined due to virtual base class
InfiniteCliqueComplex &InfiniteCliqueComplex::operator=(InfiniteCliqueComplex &&other) noexcept
{
    if (this != &other)
    {
        Network::operator=(std::move(other));
        InfiniteNetwork::operator=(std::move(other));
        CliqueComplex::operator=(std::move(other));
    }
    return *this;
}

InfiniteCliqueComplex InfiniteCliqueComplex::filter(const PointIdList &vertices) const
{
    MarkList filtered_marks{};
    filtered_marks.reserve(vertices.size());
    for (const auto &vertex : vertices)
    {
        const auto it{std::find(marks_.begin(), marks_.end(), vertex)};
        assert(it != marks_.end() && "Vertex not found");
        filtered_marks.push_back(marks_[it - marks_.begin()]);
    }
    return InfiniteCliqueComplex{
        max_dimension_,
        vertices,
        filter_edges(vertices),
        typical_vertex_mark_,
        std::move(filtered_marks)};
}

SimplexList InfiniteCliqueComplex::calc_neighbors(const Dimension dimension)
{
    // typical vertex is implicitly included in the neighbors
    if (dimension == 0U)
    {
        return SimplexList{};
    }
    if (dimension == 1U)
    {
        std::vector<Simplex> result{};
        result.reserve(vertices_.size());
        for (const auto &vertex : vertices_)
        {
            result.emplace_back(Simplex({vertex}));
        }
        return {std::move(result)};
    }
    if (dimension == 2U)
    {
        std::vector<Simplex> result{};
        result.reserve(edges_.size());
        for (const auto &edge : edges_)
        {
            result.emplace_back(Simplex({edge.first, edge.second}));
        }
        return {std::move(result)};
    }

    // implement all combinations of vertices of length dimension
    assert_simplex_tree_is_built();
    SimplexList result{};

    // iterate over simplices
    for (const auto &simplex_handle : simplex_tree_->complex_simplex_range())
    {
        if (simplex_tree_->dimension(simplex_handle) + 1 == dimension)
        {
            PointIdList vertices{};
            vertices.reserve(dimension);
            for (const auto vertex : simplex_tree_->simplex_vertex_range(simplex_handle))
            {
                vertices.push_back(vertex);
            }
            result += SimplexList{std::vector{Simplex{vertices}}};
        }
    }
    return result;
}
