#ifndef _INFINITE_NETWORK_H_
#define _INFINITE_NETWORK_H_

#include "network.h"
#include "typedefs.h"

class InfiniteNetwork : public Network
{
public:
    InfiniteNetwork(const Dimension max_dimension, const VertexId typical_vertex_id);
    SimplexHandleList get_simplices() override;
    uint32_t num_simplices() override;

private:
    VertexId typical_vertex_id_;
};

#endif