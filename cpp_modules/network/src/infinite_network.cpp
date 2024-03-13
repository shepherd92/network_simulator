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

InfiniteNetwork InfiniteNetwork::get_filtered_network(const PointIdList &vertices) const
{
    assert(std::find(vertices.begin(), vertices.end(), typical_vertex_id_) != vertices.end());
    const auto interactions{filter_simplices(interactions_, vertices)};
    InfiniteNetwork filtered_network{max_dimension_, vertices, interactions, typical_vertex_id_};
    return filtered_network;
}

Network::SimplexHandleList InfiniteNetwork::get_simplices()
{
    assert_simplicial_complex_is_built();
    const auto typical_vertex_handle{simplex_tree_->find({typical_vertex_id_})};
    return simplex_tree_->cofaces_simplex_range(typical_vertex_handle, 0U);
}

SimplexList InfiniteNetwork::get_skeleton_interactions(const Dimension max_dimension)
{
    SimplexList simplices_contatining_typical_vertex{};
    const auto simplices{get_skeleton_simplices(interactions_, max_dimension)};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            const auto vertices{simplex.vertices()};
            if (std::find(vertices.begin(), vertices.end(), typical_vertex_id_) != vertices.end())
            {
                std::lock_guard<std::mutex> lock{mutex};
                simplices_contatining_typical_vertex.push_back(simplex);
            }
        });

    sort_simplices(simplices_contatining_typical_vertex, true);
    return simplices_contatining_typical_vertex;
}

uint32_t InfiniteNetwork::num_simplices()
{
    return get_simplices().size();
}

SimplexList InfiniteNetwork::get_skeleton_simplicial_complex(const Dimension max_dimension)
{
    assert_simplicial_complex_is_built();
    const auto typical_vertex_handle{simplex_tree_->find({typical_vertex_id_})};
    const auto simplices{simplex_tree_->cofaces_simplex_range(typical_vertex_handle, max_dimension)};
    SimplexList skeleton_simplices{};
    std::mutex mutex{};

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex_handle)
        {
            const auto simplex{Simplex{get_vertices(simplex_handle)}};
            std::lock_guard<std::mutex> lock{mutex};
            skeleton_simplices.push_back(simplex);
        });

    sort_simplices(skeleton_simplices, true);
    return skeleton_simplices;
}