#ifndef _INFINITE_NETWORK_H_
#define _INFINITE_NETWORK_H_

#include "network.h"
#include "typedefs.h"

class InfiniteNetwork : public Network
{
public:
    InfiniteNetwork(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ISimplexList &interactions,
        const PointId typical_vertex_id);
    InfiniteNetwork(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const PointId typical_vertex_id);
    uint32_t num_simplices() override;
    InfiniteNetwork get_filtered_network(const PointIdList &vertices) const;

private:
    SimplexHandleList get_simplices() override;
    SimplexList get_skeleton_interactions(const Dimension max_dimension) override;
    SimplexList get_skeleton_simplicial_complex(const Dimension max_dimension) override;
    PointId typical_vertex_id_;
};

#endif