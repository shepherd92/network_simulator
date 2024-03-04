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
    const auto typical_vertex_handle{simplex_tree_->find({typical_vertex_id_})};
    return simplex_tree_->cofaces_simplex_range(typical_vertex_handle, 0U);
}

uint32_t InfiniteNetwork::num_simplices()
{
    return get_simplices().size();
}