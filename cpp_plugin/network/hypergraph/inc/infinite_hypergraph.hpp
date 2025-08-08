#ifndef _INFINITE_HYPERGRAPH_HPP_
#define _INFINITE_HYPERGRAPH_HPP_

#include "hypergraph.hpp"
#include "infinite_network.hpp"
#include "typedefs.hpp"

class InfiniteHypergraph : public InfiniteNetwork, public Hypergraph
{
public:
    InfiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const Mark typical_vertex_mark_,
        const MarkList &marks);

    InfiniteHypergraph(InfiniteHypergraph &&other) noexcept
        : Network{std::move(other)},
          InfiniteNetwork{std::move(other)},
          Hypergraph{std::move(other)}
    {
    }

    // move assignment defined due to virtual base class
    InfiniteHypergraph &operator=(InfiniteHypergraph &&other) noexcept
    {
        if (this != &other)
        {
            Network::operator=(std::move(other));
            Hypergraph::operator=(std::move(other));
            InfiniteNetwork::operator=(std::move(other));
        }
        return *this;
    }

    std::vector<Dimension> calc_interaction_dimension_distribution() const override;
    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;

private:
    SimplexList calc_neighbors(const Dimension dimension) override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;
};

#endif