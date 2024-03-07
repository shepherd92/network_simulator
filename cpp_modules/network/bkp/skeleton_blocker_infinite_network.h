#ifndef _SKELETON_BLOCKER_INFINITE_NETWORK_H_
#define _SKELETON_BLOCKER_INFINITE_NETWORK_H_

#include "network.h"
#include "skeleton_blocker_network.h"
#include "typedefs.h"

class SkeletonBlockerInfiniteNetwork : public SkeletonBlockerNetwork
{
public:
    SkeletonBlockerInfiniteNetwork(const VertexId typical_vertex_id);

    std::vector<uint32_t> calc_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;

private:
    using SkeletonBlockerSimplexRange = SkeletonBlocker::Complex_simplex_around_vertex_range;

    SkeletonBlockerSimplexRange get_simplices() const;
    VertexHandle typical_vertex_id_;
};

#endif