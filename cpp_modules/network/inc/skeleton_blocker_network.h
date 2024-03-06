#ifndef _SKELETON_BLOCKER_NETWORK_H_
#define _SKELETON_BLOCKER_NETWORK_H_

#include <gudhi/Skeleton_blocker.h>

#include "network.h"
#include "typedefs.h"

class SkeletonBlockerNetwork : public Network
{
protected:
    using SkeletonBlockerOptions = Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits;
    using SkeletonBlocker = Gudhi::skeleton_blocker::Skeleton_blocker_complex<SkeletonBlockerOptions>;

    std::optional<SkeletonBlocker> skeleton_blocker_;

private:
    virtual Complex_simplex_around_vertex_range get_simplices(const Dimension dimension) = 0;

    bool is_valid() const override;
    void reset_simplicial_complex() override;
    void initialize_simplicial_complex_if_needed() override;
    SimplexList convert_to_representable_simplices(const SimplexList &simplices_in) const override;
    void filter_simplicial_complex(const VertexList &vertices) override;
    std::vector<uint32_t> calc_degree_sequence_simplicial_complex(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;
    std::vector<Dimension> calc_simplex_dimension_distribution() override;
};

#endif // _SIMPLEX_TREE_NETWORK_H_