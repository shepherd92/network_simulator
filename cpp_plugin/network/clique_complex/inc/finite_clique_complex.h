#ifndef _FINITE_HYPERGRAPH_H_
#define _FINITE_HYPERGRAPH_H_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "finite_network.h"
#include "typedefs.h"

class FiniteCliqueComplex : public FiniteNetwork
{
public:
    FiniteCliqueComplex(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ConnectionList &edges);

private:
    void expand();
    void fill_simplicial_complex() override;
    SimplexList calc_simplices(const Dimension dimension) override;
    const ConnectionList edges_;
};

#endif