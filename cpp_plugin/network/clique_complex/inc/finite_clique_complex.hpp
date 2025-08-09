#ifndef _FINITE_CLIQUE_COMPLEX_HPP_
#define _FINITE_CLIQUE_COMPLEX_HPP_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "finite_network.hpp"
#include "typedefs.hpp"

class FiniteCliqueComplex final : public FiniteNetwork
{
public:
    FiniteCliqueComplex(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ConnectionList &edges);

    FiniteCliqueComplex(FiniteCliqueComplex &&other) noexcept;

    void set_edges(const ConnectionList &edges);
    ConnectionList get_edges() const;

    void set_vertices(const PointIdList &vertices) override;

    // move assignment defined due to virtual base class
    FiniteCliqueComplex &operator=(FiniteCliqueComplex &&other) noexcept;
    FiniteCliqueComplex filter(const PointIdList &vertices) const;

private:
    void expand();
    void fill_simplicial_complex() override;
    SimplexList calc_simplices(const Dimension dimension) override;
    ConnectionList filter_edges(const PointIdList &vertices) const;

    ConnectionList edges_;
};

#endif