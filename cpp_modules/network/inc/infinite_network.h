#ifndef _INFINITE_NETWORK_H_
#define _INFINITE_NETWORK_H_

#include "network.h"
#include "typedefs.h"

class InfiniteNetwork : public Network
{
public:
    InfiniteNetwork(
        const Dimension max_dimension,
        const VertexList &vertices,
        const ISimplexList &interactions,
        const VertexId typical_vertex_id);
    uint32_t num_simplices() override;
    InfiniteNetwork get_filtered_network(const VertexList &vertices) const;

private:
    SimplexHandleList get_simplices() override;
    SimplexList get_skeleton_interactions(const Dimension max_dimension) override;
    SimplexList get_skeleton_simplicial_complex(const Dimension max_dimension) override;
    VertexId typical_vertex_id_;
};

#endif