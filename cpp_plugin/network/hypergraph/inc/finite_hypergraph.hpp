#ifndef _FINITE_HYPERGRAPH_HPP_
#define _FINITE_HYPERGRAPH_HPP_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "hypergraph.hpp"
#include "finite_network.hpp"
#include "typedefs.hpp"

class FiniteHypergraph : public FiniteNetwork, public Hypergraph
{
public:
    FiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ISimplexList &interactions,
        const bool weighted);

    FiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const bool weighted);

    FiniteHypergraph(FiniteHypergraph &&other) noexcept
        : Network{std::move(other)},
          FiniteNetwork{std::move(other)},
          Hypergraph{std::move(other)},
          weighted_{std::move(other.weighted_)}
    {
    }

    // move assignment defined due to virtual base class
    FiniteHypergraph &operator=(FiniteHypergraph &&other) noexcept
    {
        if (this != &other)
        {
            Network::operator=(std::move(other));
            Hypergraph::operator=(std::move(other));
            FiniteNetwork::operator=(std::move(other));
            weighted_ = std::move(other.weighted_);
        }
        return *this;
    }

    std::vector<std::vector<std::pair<float, float>>> calc_persistence_intervals();
    std::vector<ISimplexList> calc_persistence_pairs();
    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;

private:
    SimplexList calc_simplices(const Dimension dimension) override;

    void fill_simplicial_complex() override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;

    bool weighted_;
};

#endif