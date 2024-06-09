#ifndef _INFINITE_NETWORK_H_
#define _INFINITE_NETWORK_H_

#include "network.h"
#include "typedefs.h"

class InfiniteNetwork : public Network
{
public:
    InfiniteNetwork(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const ISimplexList &interactions,
        const Mark typical_vertex_mark_,
        const MarkList &marks);
    InfiniteNetwork(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const Mark typical_vertex_mark_,
        const MarkList &marks);
    std::vector<Dimension> calc_facet_dimension_distribution() override;
    std::vector<Dimension> calc_interaction_dimension_distribution() const override;
    std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) override;
    std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) override;
    InfiniteNetwork filter(const PointIdList &vertices) const;
    const SimplexList &get_neighbors(const Dimension dimension);
    Mark typical_mark() const;
    SimplexList get_skeleton(const Dimension max_dimension) override;

private:
    SimplexList calc_neighbors(const Dimension dimension);
    SimplexList calc_simplices(const Dimension dimension) override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;

    Mark typical_vertex_mark_;
    MarkList marks_;
    // simplices in neighbors implicitly contain the typical vertex
    std::vector<std::optional<SimplexList>> neighbors_;
};

#endif