#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <map>
#include <optional>

#include "simplex_list.h"
#include "typedefs.h"

class Network
{
public:
    Network(const Dimension max_dimension, const PointIdList &vertices);

    Dimension get_max_dimension() const;
    void set_max_dimension(const Dimension dimension);
    uint32_t num_simplices(const Dimension dimension);

    virtual void reset();

    std::vector<Dimension> calc_simplex_dimension_distribution();
    virtual std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) = 0;

    virtual void keep_only_vertices(const PointIdList &vertices);
    ISimplexList get_skeleton_interface(const Dimension max_dimension);

    uint32_t num_vertices();
    void set_vertices(const PointIdList &vertices);
    PointIdList get_vertices() const;

protected:
    virtual SimplexList calc_simplices(const Dimension dimension) = 0;
    virtual SimplexList get_skeleton(const Dimension max_dimension) = 0;
    const SimplexList &get_simplices(const Dimension dimension);

    Dimension max_dimension_;
    std::vector<std::optional<SimplexList>> simplices_;

private:
    PointIdList vertices_;
    std::map<PointId, std::vector<PointId>> create_vertex_simplex_map(const SimplexList &simplices) const;
};

#endif
