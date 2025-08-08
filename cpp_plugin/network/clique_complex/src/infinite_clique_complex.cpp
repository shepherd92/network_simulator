#include <cassert>
#include <mutex>

#include "infinite_clique_complex.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"

InfiniteCliqueComplex::InfiniteCliqueComplex(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices},
      InfiniteNetwork{typical_vertex_mark_, marks}
{
}

InfiniteCliqueComplex::InfiniteCliqueComplex(InfiniteCliqueComplex &&other) noexcept
    : Network{std::move(other)},
      InfiniteNetwork{std::move(other)}
{
}

// move assignment defined due to virtual base class
InfiniteCliqueComplex &InfiniteCliqueComplex::operator=(InfiniteCliqueComplex &&other) noexcept
{
    if (this != &other)
    {
        Network::operator=(std::move(other));
        InfiniteNetwork::operator=(std::move(other));
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
    return InfiniteCliqueComplex{max_dimension_, vertices, typical_vertex_mark_, std::move(filtered_marks)};
}

SimplexList InfiniteCliqueComplex::calc_neighbors(const Dimension dimension)
{
    if (dimension == 0U)
    {
        return SimplexList{};
    }
    // typical vertex is implicitly included in the neighbors
    // implement all combinations of vertices of length dimension
    Simplex simplex_all_vertices{vertices_};
    const auto neighbors{simplex_all_vertices.faces(dimension)};
    return neighbors;
}
