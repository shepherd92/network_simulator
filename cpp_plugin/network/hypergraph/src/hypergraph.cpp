#include <atomic>
#include <cassert>
#include <execution>
#include <mutex>

#include "hypergraph.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"
#include "typedefs.hpp"

// Interface functions exposed to Python

void Hypergraph::set_interactions(const ISimplexList &interactions)
{
    interactions_ = SimplexList{interactions};
}

ISimplexList Hypergraph::get_interactions() const
{
    return interactions_.raw();
}

// Internal functions

Hypergraph::Hypergraph(const SimplexList &interactions)
    : interactions_{interactions}
{
}

Hypergraph::Hypergraph(Hypergraph &&other) noexcept
{
    interactions_ = std::move(other.interactions_);
}

Hypergraph &Hypergraph::operator=(Hypergraph &&other) noexcept
{
    if (this != &other)
    {
        interactions_ = std::move(other.interactions_);
    }
    return *this;
}

std::vector<Dimension> Hypergraph::calc_interaction_dimension_distribution() const
{
    return interactions_.calc_dimension_distribution();
}

void Hypergraph::keep_only_vertices(const PointIdList &vertices)
{
    Network::keep_only_vertices(vertices);
    interactions_ = interactions_.filter(vertices);
}
