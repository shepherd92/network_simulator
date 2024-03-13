#ifndef _INFINITE_HYPERGRAPH_MODEL_H_
#define _INFINITE_HYPERGRAPH_MODEL_H_

#include "hypergraph_model.h"

class InfiniteNetwork;

class InfiniteHypergraphModel : public HypergraphModel
{
public:
    InfiniteHypergraphModel(const py::array_t<double> &parameters_in, const uint32_t seed);
    std::vector<InfiniteNetwork> generate_networks(const uint32_t num_of_networks) const;
    InfiniteNetwork generate_network() const;

private:
    NetworkInterface generate_finite_network_interface(
        const py::array_t<double> &model_parameters_input,
        const uint32_t seed);

    PointList create_interactions(const Mark u) const;
    PointList create_points(const size_t num_of_nodes, const float exponent) const;
    PointList create_vertices(const PointList &interactions) const;
    PositionList generate_positions(
        const MarkList &marks,
        const float beta_x_u_to_gamma,
        const float exponent) const;
    PointList create_vertices_in_interaction_neighborhood(const Point &interaction) const;
    PositionList generate_positions_in_interaction_neighborhood(const Point &interaction, const MarkList &marks) const;
    PositionList generate_positions_in_vertex_neighborhood(const Point &vertex, const MarkList &marks) const;
    PositionList generate_positions_in_neighborhood(
        const Point &point,
        const MarkList &marks,
        const float exponent_of_central_point,
        const float exponent_of_points_in_neighborhood) const;

    float distance(const Point &first, const Point &second) const override;
};

#endif