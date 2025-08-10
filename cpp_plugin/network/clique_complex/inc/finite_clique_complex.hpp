#ifndef _FINITE_CLIQUE_COMPLEX_HPP_
#define _FINITE_CLIQUE_COMPLEX_HPP_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "clique_complex.hpp"
#include "finite_network.hpp"
#include "typedefs.hpp"

class FiniteCliqueComplex final : public FiniteNetwork, public CliqueComplex
{
public:
    FiniteCliqueComplex(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ConnectionList &edges);

    FiniteCliqueComplex(FiniteCliqueComplex &&other) noexcept;

    // move assignment defined due to virtual base class
    FiniteCliqueComplex &operator=(FiniteCliqueComplex &&other) noexcept;
    FiniteCliqueComplex filter(const PointIdList &vertices) const;

private:
    SimplexList calc_simplices(const Dimension dimension) override;
};

#endif