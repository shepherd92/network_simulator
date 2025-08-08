#ifndef _INFINITE_CLIQUE_COMPLEX_H_
#define _INFINITE_CLIQUE_COMPLEX_H_

#include "infinite_network.h"
#include "typedefs.h"

class InfiniteCliqueComplex : public InfiniteNetwork
{
public:
    InfiniteCliqueComplex(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const Mark typical_vertex_mark_,
        const MarkList &marks);

private:
    SimplexList calc_neighbors(const Dimension dimension) override;
};

#endif