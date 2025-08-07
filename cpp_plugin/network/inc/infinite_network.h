#ifndef _INFINITE_NETWORK_H_
#define _INFINITE_NETWORK_H_

#include "network.h"
#include "typedefs.h"

class InfiniteNetwork : virtual public Network
{
public:
    InfiniteNetwork(
        const Mark typical_vertex_mark_,
        const MarkList &marks);
    std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;
    const SimplexList &get_neighbors(const Dimension dimension);
    Mark typical_mark() const;
    SimplexList get_skeleton(const Dimension max_dimension) override;

private:
    SimplexList calc_simplices(const Dimension dimension) override;
    virtual SimplexList calc_neighbors(const Dimension dimension) = 0;

    const Mark typical_vertex_mark_;
    MarkList marks_;
    // simplices in neighbors implicitly contain the typical vertex
    std::vector<std::optional<SimplexList>> neighbors_;
};

#endif