#ifndef _HYPERGRAPH_H_
#define _HYPERGRAPH_H_

#include <map>
#include <optional>

#include "network.h"
#include "simplex_list.h"
#include "typedefs.h"

class Hypergraph : virtual public Network
{
public:
    Hypergraph(const SimplexList &interactions);
    virtual std::vector<Dimension> calc_interaction_dimension_distribution() const;
    virtual std::vector<uint32_t> calc_simplex_interaction_degree_sequence(
        const Dimension simplex_dimension) = 0;

    ISimplexList get_interactions_interface() const;
    void set_interactions(const ISimplexList &interactions);

    void keep_only_vertices(const PointIdList &vertices) override;

protected:
    const SimplexList &get_interactions() const;

private:
    SimplexList interactions_;
    virtual std::vector<uint32_t> calc_vertex_interaction_degree_distribution() const = 0;
};

#endif
