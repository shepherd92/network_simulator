#include <mutex>

#include "infinite_network.h"

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ISimplexList &interactions,
    const PointId typical_vertex_id)
    : Network{max_dimension, vertices, create_simplices(interactions)},
      typical_vertex_id_{typical_vertex_id}
{
}

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &interactions,
    const PointId typical_vertex_id)
    : Network{max_dimension, vertices, interactions},
      typical_vertex_id_{typical_vertex_id}
{
}

PointIdList InfiniteNetwork::get_vertices() const
{
    return PointIdList{typical_vertex_id_};
}

SimplexList InfiniteNetwork::calc_simplices(const Dimension dimension)
{
    const auto all_simplices_of_dimension{get_faces_simplices(get_facets(), dimension)};
    return get_cofaces_simplices(all_simplices_of_dimension, Simplex{PointIdList{typical_vertex_id_}});
}

InfiniteNetwork InfiniteNetwork::filter(const PointIdList &vertices) const
{
    InfiniteNetwork result{*this};
    result.keep_only_vertices(vertices);
    return result;
}