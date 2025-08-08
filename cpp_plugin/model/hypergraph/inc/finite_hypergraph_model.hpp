#ifndef _FINITE_HYPERGRAPH_MODEL_HPP_
#define _FINITE_HYPERGRAPH_MODEL_HPP_

#include "finite_model.hpp"
#include "hypergraph_model.hpp"

class FiniteHypergraph;

class FiniteHypergraphModel : public FiniteModel, public HypergraphModel
{
public:
    FiniteHypergraphModel(const std::vector<double> &parameters_in, const uint32_t seed);
    std::tuple<FiniteHypergraph, MarkPositionList, MarkPositionList> generate_network() const;

private:
    bool rectangle_points_surely_connect(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const override;
    const bool weighted_;
};

#endif