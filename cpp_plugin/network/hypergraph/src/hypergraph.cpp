#include <atomic>
#include <cassert>
#include <execution>
#include <mutex>

#include "hypergraph.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"
#include "typedefs.h"

Hypergraph::Hypergraph(const SimplexList &interactions)
    : interactions_{interactions}
{
}

ISimplexList Hypergraph::get_interactions_interface() const
{
    return get_interactions().raw();
}

std::vector<Dimension> Hypergraph::calc_interaction_dimension_distribution() const
{
    return get_interactions().calc_dimension_distribution();
}

const SimplexList &Hypergraph::get_interactions() const
{
    return interactions_;
}

void Hypergraph::set_interactions(const SimplexList &interactions)
{
    interactions_ = SimplexList{interactions};
}

void Hypergraph::keep_only_vertices(const PointIdList &vertices)
{
    Network::keep_only_vertices(vertices);
    set_interactions(get_interactions().filter(vertices));
}
