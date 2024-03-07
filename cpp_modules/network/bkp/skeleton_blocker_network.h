#ifndef _SKELETON_BLOCKER_NETWORK_H_
#define _SKELETON_BLOCKER_NETWORK_H_

#include <gudhi/Skeleton_blocker.h>
#include <optional>

#include "network.h"
#include "typedefs.h"

class SkeletonBlockerNetwork : public Network
{
public:
protected:
    using SkeletonBlockerOptions = Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits;
    using SkeletonBlocker = Gudhi::skeleton_blocker::Skeleton_blocker_complex<SkeletonBlockerOptions>;
    using SkeletonBlockerSimplex = SkeletonBlocker::Simplex;
    using VertexHandle = SkeletonBlocker::Vertex_handle;

    bool is_valid() const override;

    std::optional<SkeletonBlocker> skeleton_blocker_;

private:
    template <size_t num_args>
    struct create_skeleton_blocker_simplex
    {
    public:
        void operator()(const VertexList &vertices) const
        {
            call(vertices, BuildIndices<num_args>{});
        }

    private:
        template <typename size_t... I>
        void call(const VertexList &vertices, indices<I...>)
        {
            SkeletonBlockerSimplex(vertices[I]...);
        }
    };

    void add_simplex(const VertexList &simplex) override;
    void reset_simplicial_complex() override;
    void initialize_simplicial_complex_if_needed() override;
    SimplexList convert_to_representable_simplices(const SimplexList &simplices_in) const override;
    void filter_simplicial_complex(const VertexList &vertices) override;
};

#endif // _SIMPLEX_TREE_NETWORK_H_