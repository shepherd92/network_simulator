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
    const PointId typical_vertex_id,
    const MarkList &marks)
    : Network{max_dimension, vertices, SimplexList{interactions}},
      typical_vertex_id_{typical_vertex_id},
      marks_{marks},
      neighbors_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
}

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &interactions,
    const PointId typical_vertex_id,
    const MarkList &marks)
    : Network{max_dimension, vertices, interactions},
      typical_vertex_id_{typical_vertex_id},
      marks_{marks},
      neighbors_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
}

PointIdList InfiniteNetwork::get_vertices() const
{
    return PointIdList{typical_vertex_id_};
}

SimplexList InfiniteNetwork::calc_simplices(const Dimension dimension)
{
    const auto cofaces_of_typical_vertex{get_neighbors(dimension)};

    const auto typical_vertex_mark{marks_[typical_vertex_id_]};
    std::vector<Simplex> result{};
    for (auto &simplex : cofaces_of_typical_vertex)
    {
        bool typical_vertex_is_oldest{true};
        for (const auto &vertex : simplex.vertices())
        {
            if (vertex != typical_vertex_id_ && marks_[vertex] < typical_vertex_mark)
            {
                typical_vertex_is_oldest = false;
                break;
            }
        }
        if (typical_vertex_is_oldest)
        {
            result.push_back(simplex);
        }
    }
    return SimplexList{result};
}

SimplexList InfiniteNetwork::calc_neighbors(const Dimension dimension)
{
    const auto all_simplices_of_dimension{get_facets().faces(dimension)};
    return all_simplices_of_dimension.cofaces(Simplex{PointIdList{typical_vertex_id_}});
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

InfiniteNetwork InfiniteNetwork::filter(const PointIdList &vertices) const
{
    InfiniteNetwork result{*this};
    result.keep_only_vertices(vertices);
    return result;
}