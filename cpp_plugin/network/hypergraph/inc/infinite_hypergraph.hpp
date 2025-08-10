#ifndef _INFINITE_HYPERGRAPH_HPP_
#define _INFINITE_HYPERGRAPH_HPP_

#include "hypergraph.hpp"
#include "infinite_network.hpp"
#include "typedefs.hpp"

class InfiniteHypergraph final : public InfiniteNetwork, public Hypergraph
{
public:
    InfiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const Mark typical_vertex_mark_,
        const MarkList &marks);

    InfiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ISimplexList &interactions,
        const Mark typical_vertex_mark_,
        const MarkList &marks);

    // move defined due to virtual base class
    InfiniteHypergraph(InfiniteHypergraph &&other) noexcept;
    InfiniteHypergraph &operator=(InfiniteHypergraph &&other) noexcept;

    std::vector<Dimension> calc_interaction_dimension_distribution() const override;
    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;

    InfiniteHypergraph filter(const PointIdList &vertices);

    void fill_simplicial_complex() override;

private:
    SimplexList calc_neighbors(const Dimension dimension) override;
};

#endif