#ifndef _INFINITE_CLIQUE_COMPLEX_HPP_
#define _INFINITE_CLIQUE_COMPLEX_HPP_

#include "infinite_network.hpp"
#include "typedefs.hpp"

class InfiniteCliqueComplex final : public InfiniteNetwork
{
public:
    InfiniteCliqueComplex(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const Mark typical_vertex_mark_,
        const MarkList &marks);

    // move assignment defined due to virtual base class
    InfiniteCliqueComplex(InfiniteCliqueComplex &&other) noexcept;
    InfiniteCliqueComplex &operator=(InfiniteCliqueComplex &&other) noexcept;

    InfiniteCliqueComplex filter(const PointIdList &vertices) const;

private:
    SimplexList calc_neighbors(const Dimension dimension) override;
};

#endif