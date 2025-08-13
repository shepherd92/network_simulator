#ifndef _HYPERGRAPH_MODEL_INL_
#define _HYPERGRAPH_MODEL_INL_

#include "point.hpp"
#include "typedefs.hpp"

bool HypergraphModel::connects(const Point &vertex, const Point &interaction) const
{
    return distance(vertex, interaction) < vertex.mark() * interaction.mark();
}

Dimension HypergraphModel::max_dimension() const
{
    return parameters_.max_dimension;
}

float HypergraphModel::beta() const
{
    return parameters_.beta;
}

float HypergraphModel::gamma() const
{
    return parameters_.gamma;
}

float HypergraphModel::gamma_prime() const
{
    return parameters_.gamma_prime;
}

float HypergraphModel::lambda() const
{
    return parameters_.network_size;
}

float HypergraphModel::lambda_prime() const
{
    return parameters_.interaction_intensity;
}

#endif