#ifndef _FINITE_HYPERGRAPH_HPP_
#define _FINITE_HYPERGRAPH_HPP_

#include <iostream>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include "hypergraph.hpp"
#include "finite_network.hpp"
#include "typedefs.hpp"

class FiniteHypergraph final : public FiniteNetwork, public Hypergraph
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

    // move defined due to virtual base class
    FiniteHypergraph(FiniteHypergraph &&other) noexcept;
    FiniteHypergraph &operator=(FiniteHypergraph &&other) noexcept;

    std::vector<std::vector<std::pair<float, float>>> calc_persistence_intervals();
    std::vector<ISimplexList> calc_persistence_pairs();
    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;

    FiniteHypergraph filter(const PointIdList &vertices);

private:
    SimplexList calc_simplices(const Dimension dimension) override;

    void fill_simplicial_complex() override;

    bool weighted_;
};

#endif