#ifndef _FINITE_HYPERGRAPH_MODEL_H_
#define _FINITE_HYPERGRAPH_MODEL_H_

#include "hypergraph_model.h"

class FiniteNetwork;

class FiniteHypergraphModel : public HypergraphModel
{
public:
    FiniteHypergraphModel(const py::array_t<double> &parameters_in, const uint32_t seed);
    FiniteNetwork generate_network() const;

private:
    NetworkInterface generate_finite_network_interface(
        const py::array_t<double> &model_parameters_input,
        const uint32_t seed);

    PointList create_points(const size_t num_of_nodes) const;
    PositionList generate_positions(const size_t num_of_nodes) const;
    float distance(const Point &first, const Point &second) const override;
};

#endif