#ifndef _INFINITE_HYPERGRAPH_MODEL_H_
#define _INFINITE_HYPERGRAPH_MODEL_H_

#include "hypergraph_model.h"
#include "infinite_model.h"

class Center;
class Hyperbola;
class InfiniteNetwork;
class NeighborhoodPart;

class InfiniteHypergraphModel : public InfiniteModel, public HypergraphModel
{
public:
    InfiniteHypergraphModel(const std::vector<double> &parameters_in, const uint32_t seed);

    std::vector<std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList>>
    generate_networks(const uint32_t num_of_infinite_networks) const;

private:
    std::tuple<InfiniteNetwork, MarkPositionList, MarkPositionList> generate_network() const;
    PointList create_interactions(const Mark u) const;
    PointList transform_interactions(const PointList &interactions) const;

    std::pair<std::vector<Center>, std::vector<Hyperbola>>
    create_dominating_neighborhood_parts(const PointList &transformed_interactions) const;

    std::vector<Center> create_neighborhood_centers(const PointList &transformed_interactions) const;
    std::vector<Center> merge_neighborhood_centers(std::vector<Center> &centers) const;
    std::vector<Hyperbola> create_neighborhood_tails(const PointList &transformed_interactions) const;
    std::vector<std::vector<Hyperbola>> remove_center_domains_from_neighborhood_tails(
        const std::vector<Hyperbola> &tails,
        const std::vector<Center> &centers) const;
    std::vector<Hyperbola> determine_dominating_hyperbolas(
        const std::vector<std::vector<Hyperbola>> &tails_in_slots) const;
    PointList create_vertices(
        const std::vector<Center> &centers,
        const std::vector<Hyperbola> &hyperbolas) const;
    PointList create_vertices_under_centers(const std::vector<Center> &centers) const;
    PointList create_vertices_under_tails(const std::vector<Hyperbola> &hyperbolas) const;

    Hyperbola get_neighborhood_left_tail(const Point &point) const;
    Center get_neighborhood_center(const Point &point) const;
    Hyperbola get_neighborhood_right_tail(const Point &point) const;

    PointList create_vertices_in_interaction_neighborhood(const Point &interaction) const;
    PositionList generate_positions_in_vertex_neighborhood(const Point &vertex, const MarkList &marks) const;

    bool rectangle_points_surely_connect(
        const Rectangle &vertex_rectangle,
        const Rectangle &interaction_rectangle) const override;
};

#endif