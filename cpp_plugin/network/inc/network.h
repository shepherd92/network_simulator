#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <map>
#include <optional>

#include "simplex_list.h"
#include "typedefs.h"

class Network
{
public:
    Network(
        const Dimension max_dimension,
        const PointIdList &vertices,
        const SimplexList &interactions);

    Dimension get_max_dimension() const;
    void set_max_dimension(const Dimension dimension);
    uint32_t num_simplices(const Dimension dimension);

    void keep_only_vertices(const PointIdList &vertices);
    virtual void reset();

    std::vector<Dimension> calc_simplex_dimension_distribution();
    virtual std::vector<Dimension> calc_facet_dimension_distribution();
    virtual std::vector<Dimension> calc_interaction_dimension_distribution() const;
    virtual std::vector<uint32_t> calc_coface_degree_sequence(
        const Dimension simplex_dimension,
        const Dimension neighbor_dimension) = 0;
    virtual std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) = 0;

    ISimplexList get_interactions_interface() const;
    ISimplexList get_facets_interface();
    ISimplexList get_skeleton_interface(const Dimension max_dimension);

    uint32_t num_vertices();
    void set_vertices(const PointIdList &vertices);

    void set_interactions(const ISimplexList &interactions);

protected:
    virtual SimplexList get_skeleton(const Dimension max_dimension) = 0;

    const SimplexList &get_interactions() const;
    const SimplexList &get_facets();
    const SimplexList &get_simplices(const Dimension dimension);

    Dimension max_dimension_;
    PointIdList vertices_;
    SimplexList interactions_;
    std::optional<SimplexList> facets_;
    std::vector<std::optional<SimplexList>> simplices_;

private:
    virtual SimplexList calc_simplices(const Dimension dimension) = 0;
    virtual std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const = 0;
    std::map<PointId, std::vector<int32_t>> create_vertex_simplex_map(const SimplexList &simplices) const;
    void filter_facets(const PointIdList &vertices);
};

#endif