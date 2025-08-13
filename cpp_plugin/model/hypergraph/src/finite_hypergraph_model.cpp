#include <random>

#include "finite_hypergraph_model.hpp"
#include "finite_hypergraph.hpp"
#include "point.hpp"
#include "simplex_list.hpp"
#include "rectangle.hpp"

FiniteHypergraphModel::FiniteHypergraphModel(const std::vector<double> &parameters_in, const uint32_t seed)
    : Model{seed}, FiniteModel{seed}, HypergraphModel{parameters_in}, weighted_{parameters_in[6] > 0.5}
{
}

std::tuple<FiniteHypergraph, MarkPositionList, MarkPositionList> FiniteHypergraphModel::generate_network() const
{
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(lambda() * torus_size)(random_number_generator_)};
    const auto vertices{create_points(num_of_vertices)};
    const auto vertex_ids{convert_to_id_list(vertices)};
    auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    const auto num_of_interactions{std::poisson_distribution<uint32_t>(lambda_prime() * torus_size)(random_number_generator_)};
    const auto interactions{create_points(num_of_interactions)};
    const auto interaction_ids{convert_to_id_list(interactions)};
    auto interaction_mark_position_pairs{convert_to_mark_position_pairs(interactions)};

    const auto connections{generate_connections(vertices, interactions)};
    const auto interaction_simplices{create_interaction_simplices_from_connections(interaction_ids, connections)};
    FiniteHypergraph network{max_dimension(), vertex_ids, interaction_simplices, weighted_};

    return std::tuple<FiniteHypergraph, MarkPositionList, MarkPositionList>(
        std::move(network),
        std::move(vertex_mark_position_pairs),
        std::move(interaction_mark_position_pairs));
}

bool FiniteHypergraphModel::rectangle_points_surely_connect(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const
{
    const Point vtl{vertex_rectangle.top(), vertex_rectangle.left()};
    const Point vtr{vertex_rectangle.top(), vertex_rectangle.right()};
    const Point itl{interaction_rectangle.top(), interaction_rectangle.left()};
    const Point itr{interaction_rectangle.top(), interaction_rectangle.right()};

    // all points connect within the torus (no wrapping around the torus boundaries)
    if (fabs(vtl.position() - itl.position()) < vtl.mark() * itl.mark() &&
        fabs(vtl.position() - itr.position()) < vtl.mark() * itr.mark() &&
        fabs(vtr.position() - itl.position()) < vtr.mark() * itl.mark() &&
        fabs(vtr.position() - itr.position()) < vtr.mark() * itr.mark())
    {
        return true;
    }

    // all points connect wrapping around the torus
    if (torus_size - fabs(vtl.position() - itl.position()) < vtl.mark() * itl.mark() &&
        torus_size - fabs(vtl.position() - itr.position()) < vtl.mark() * itr.mark() &&
        torus_size - fabs(vtr.position() - itl.position()) < vtr.mark() * itl.mark() &&
        torus_size - fabs(vtr.position() - itr.position()) < vtr.mark() * itr.mark())
    {
        return true;
    }

    // assumption: no rectangles wrap around the torus boundaries
    return false;
}
