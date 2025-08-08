#ifndef _FINITE_HYPERGRAPH_HPP_
#define _FINITE_HYPERGRAPH_HPP_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "clique_complex.hpp"
#include "finite_network.hpp"
#include "typedefs.hpp"

class FiniteCliqueComplex : public FiniteNetwork, CliqueComplex
{
public:
    FiniteCliqueComplex(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ConnectionList &edges);

    FiniteCliqueComplex(FiniteCliqueComplex &&other) noexcept
        : Network{std::move(other)},
          FiniteNetwork{std::move(other)},
          edges_{std::move(other.edges_)}
    {
    }

    // move assignment defined due to virtual base class
    FiniteCliqueComplex &operator=(FiniteCliqueComplex &&other) noexcept
    {
        if (this != &other)
        {
            Network::operator=(std::move(other));
            FiniteNetwork::operator=(std::move(other));
            edges_ = std::move(other.edges_);
        }
        return *this;
    }

private:
    void expand();
    void fill_simplicial_complex() override;
    SimplexList calc_simplices(const Dimension dimension) override;
    ConnectionList edges_;
};

#endif