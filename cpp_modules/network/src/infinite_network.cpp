#include <mutex>

#include "infinite_network.h"

InfiniteNetwork::InfiniteNetwork(
    const VertexList &vertices,
    const ISimplexList &interactions,
    const VertexId typical_vertex_id)
    : Network{vertices, interactions},
      typical_vertex_id_{typical_vertex_id}
{
}

Network::SimplexHandleList InfiniteNetwork::get_simplices()
{
    assert(is_valid());
    const auto typical_vertex_handle{simplex_tree_->find({typical_vertex_id_})};
    return simplex_tree_->cofaces_simplex_range(typical_vertex_handle, 0U);
}

uint32_t InfiniteNetwork::num_simplices()
{
    return get_simplices().size();
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
    return simplices_contatining_typical_vertex;
}

SimplexList InfiniteNetwork::get_skeleton_simplicial_complex(const Dimension max_dimension)
{
    assert(is_valid());
    const auto typical_vertex_handle{simplex_tree_->find({typical_vertex_id_})};
    auto simplices{simplex_tree_->cofaces_simplex_range(typical_vertex_handle, max_dimension)};
    SimplexList skeleton_simplices{};

    std::transform(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        std::back_inserter(skeleton_simplices),
        [this](const auto &simplex_handle) -> Simplex
        {
            return Simplex{get_vertices(simplex_handle)};
        });
    return skeleton_simplices;
}