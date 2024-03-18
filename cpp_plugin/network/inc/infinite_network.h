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
    PointIdList get_vertices() const override;
    InfiniteNetwork filter(const PointIdList &vertices) const;

private:
    SimplexList calc_simplices(const Dimension dimension) override;
    PointId typical_vertex_id_;
};

#endif