#ifndef _INFINITE_HYPERGRAPH_H_
#define _INFINITE_HYPERGRAPH_H_

#include "hypergraph.h"
#include "infinite_network.h"
#include "typedefs.h"

class InfiniteHypergraph : public InfiniteNetwork, public Hypergraph
{
public:
    InfiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ISimplexList &interactions,
        const Mark typical_vertex_mark_,
        const MarkList &marks);
    InfiniteHypergraph(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const Mark typical_vertex_mark_,
        const MarkList &marks);
    std::vector<Dimension> calc_interaction_dimension_distribution() const override;
    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;
    const SimplexList &get_neighbors(const Dimension dimension);

private:
    SimplexList calc_neighbors(const Dimension dimension) override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;
};

#endif