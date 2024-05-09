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
        const PointId typical_vertex_id,
        const MarkList &marks);
    InfiniteNetwork(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions,
        const PointId typical_vertex_id,
        const MarkList &marks);
    PointIdList get_vertices() const override;
    std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const override;
    InfiniteNetwork filter(const PointIdList &vertices) const;
    const SimplexList &get_neighbors(const Dimension dimension) override;

private:
    SimplexList calc_neighbors(const Dimension dimension);
    SimplexList calc_simplices(const Dimension dimension) override;
    Simplex get_typical_vertex_as_simplex() const;

    PointId typical_vertex_id_;
    MarkList marks_;
    std::vector<std::optional<SimplexList>> neighbors_;
};

#endif