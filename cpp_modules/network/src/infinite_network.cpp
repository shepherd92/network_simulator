#include "infinite_network.h"

InfiniteNetwork::InfiniteNetwork(
    const Dimension max_dimension,
    const VertexId typical_vertex_id)
    : Network{max_dimension},
      typical_vertex_id_{typical_vertex_id}
{
}

SimplexHandleList InfiniteNetwork::get_simplices()
{
    switch (type_)
    {
    case SimplicialComplexType::simplex_tree:
        assert(simplex_tree_.has_value());
        const auto typical_vertex_handle{simplex_tree_->find({typical_vertex_id_})};
        return simplex_tree_->cofaces_simplex_range(typical_vertex_handle, 0U);
    case SimplicialComplexType::skeleton_blocker:
        assert(skeleton_blocker_.has_value());
        skeleton_blocker_
            skeleton_blocker_->;
        return skeleton_blocker_->star_simplex_range(typical_vertex_id_);
    default:
        assert(false);
    }
}

uint32_t InfiniteNetwork::num_simplices()
{
    return get_simplices().size();
}