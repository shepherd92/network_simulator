#ifndef _SIMPLEX_TREE_INFINITE_NETWORK_H_
#define _SIMPLEX_TREE_INFINITE_NETWORK_H_

#include "simplex_tree_network.h"
#include "typedefs.h"

class SimplexTreeInfiniteNetwork : public SimplexTreeNetwork
{
public:
    uint32_t num_simplices() override;

private:
    SimplexHandleList get_simplices() override;
    VertexId typical_vertex_id_;
};

uint32_t SimplexTreeInfiniteNetwork::num_simplices()
{
    return get_simplices().size();
}

SimplexTreeNetwork::SimplexHandleList SimplexTreeInfiniteNetwork::get_simplices()
{
    assert(is_valid());
    const auto typical_vertex_handle{simplex_tree_->find({typical_vertex_id_})};
    return simplex_tree_->cofaces_simplex_range(typical_vertex_handle, 0U);
}

#endif