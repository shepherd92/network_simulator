#ifndef _INFINITE_HYPERGRAPH_MODEL_H_
#define _INFINITE_HYPERGRAPH_MODEL_H_

#include "hypergraph_model.h"
#include "infinite_model.h"

class Center;
class Hyperbola;
class InfiniteNetwork;

class InfiniteHypergraphModel : public InfiniteModel, public HypergraphModel
{
public:
    InfiniteHypergraphModel(const std::vector<double> &parameters_in, const uint32_t seed);

    std::vector<std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>>
    generate_networks(const uint32_t num_of_infinite_networks) const;

    Hyperbola get_neighborhood_left_tail(const Point &point) const;
    Center get_neighborhood_center(const Point &point) const;
    Hyperbola get_neighborhood_right_tail(const Point &point) const;

private:
    std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList> generate_network() const;
    PointList create_interactions(const Mark u) const;
    PointList transform_interactions(const PointList &interactions) const;

    std::vector<std::unique_ptr<NeighborhoodPart>>
    create_dominating_neighborhood_parts(const PointList &transformed_interactions) const;

    std::vector<Center> create_neighborhood_centers(const PointList &transformed_interactions) const;
    void merge_neighborhood_centers(std::vector<Center> &centers) const;
    std::vector<Hyperbola> create_neighborhood_tails(const PointList &transformed_interactions) const;
    void remove_center_domains_from_neighborhood_tails(
        std::vector<Hyperbola> &tails,
        const std::vector<Center> &centers) const;
    std::vector<Hyperbola> determine_dominating_hyperbolas(std::vector<Hyperbola> &tails) const;

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

    bool rectangle_points_surely_connect(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const override;
};

#endif