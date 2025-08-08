#ifndef _INFINITE_CLIQUE_COMPLEX_HPP_
#define _INFINITE_CLIQUE_COMPLEX_HPP_

#include "infinite_network.hpp"
#include "typedefs.hpp"

class InfiniteCliqueComplex : public InfiniteNetwork
{
public:
    InfiniteCliqueComplex(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const Mark typical_vertex_mark_,
        const MarkList &marks);

    InfiniteCliqueComplex(InfiniteCliqueComplex &&other) noexcept
        : Network{std::move(other)},
          InfiniteNetwork{std::move(other)}
    {
    }

    // move assignment defined due to virtual base class
    InfiniteCliqueComplex &operator=(InfiniteCliqueComplex &&other) noexcept
    {
        if (this != &other)
        {
            Network::operator=(std::move(other));
            InfiniteNetwork::operator=(std::move(other));
        }
        return *this;
    }

private:
    SimplexList calc_neighbors(const Dimension dimension) override;
};

#endif