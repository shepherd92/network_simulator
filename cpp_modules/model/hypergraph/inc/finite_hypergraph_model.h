#ifndef _FINITE_HYPERGRAPH_MODEL_H_
#define _FINITE_HYPERGRAPH_MODEL_H_

#include "finite_model.h"
#include "hypergraph_model.h"

class FiniteNetwork;

class FiniteHypergraphModel : public FiniteModel, public HypergraphModel
{
public:
    FiniteHypergraphModel(const py::array_t<double> &parameters_in, const uint32_t seed);
    std::tuple<FiniteNetwork, MarkPositionList, MarkPositionList> generate_network() const;
};

#endif