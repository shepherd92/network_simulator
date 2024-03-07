#ifndef _SKELETON_BLOCKER_FINITE_NETWORK_H_
#define _SKELETON_BLOCKER_FINITE_NETWORK_H_

#include "skeleton_blocker_network.h"
#include "typedefs.h"

class SkeletonBlockerFiniteNetwork : public SkeletonBlockerNetwork
{
public:
    std::vector<uint32_t> calc_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;

private:
    template <typename SimplexIterator>
    std::vector<uint32_t> calc_degree_sequence(
        const SimplexIterator &simplices,
        const Dimension neighbor_dimension);
};

#endif