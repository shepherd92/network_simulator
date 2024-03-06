#ifndef _SKELETON_BLOCKER_FINITE_NETWORK_H_
#define _SKELETON_BLOCKER_FINITE_NETWORK_H_

#include "skeleton_blocker_network.h"
#include "typedefs.h"

class SkeletonBlockerFiniteNetwork : public SkeletonBlockerNetwork
{
public:
    SkeletonBlockerFiniteNetwork();
    SimplexHandleList get_simplices(const Dimension dimension) override;
    uint32_t num_simplices() override;
};

#endif