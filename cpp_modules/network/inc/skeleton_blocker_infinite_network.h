#ifndef _SKELETON_BLOCKER_INFINITE_NETWORK_H_
#define _SKELETON_BLOCKER_INFINITE_NETWORK_H_

#include "network.h"
#include "skeleton_blocker_network.h"
#include "typedefs.h"

class SkeletonBlockerInfiniteNetwork : public SkeletonBlockerNetwork
{
public:
    SkeletonBlockerInfiniteNetwork();

private:
    Complex_simplex_around_vertex_range get_simplices() override;
    VertexId typical_vertex_id_;
};

SkeletonBlockerNetwork::SimplexHandleList SkeletonBlockerInfiniteNetwork::get_simplices()
{
    assert(is_valid());
    return skeleton_blocker_->star_simplex_range(typical_vertex_id_);
}

#endif