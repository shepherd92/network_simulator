#include <random>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "finite_hypergraph_model.h"
#include "finite_network.h"
#include "point.h"
#include "rectangle.h"

namespace py = pybind11;

constexpr float TORUS_SIZE{1.};

FiniteHypergraphModel::FiniteHypergraphModel(const py::array_t<double> &parameters_in, const uint32_t seed)
    : HypergraphModel{parameters_in, seed}
{
}

FiniteNetwork FiniteHypergraphModel::generate_network() const
{
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(lambda() * TORUS_SIZE)(random_number_generator_)};
    const auto vertices{create_points(num_of_vertices)};
    const auto vertex_ids{convert_to_id_list(vertices)};
    const auto vertex_mark_position_pairs{convert_to_mark_position_pairs(vertices)};

    const auto num_of_interactions{std::poisson_distribution<uint32_t>(lambda_prime() * TORUS_SIZE)(random_number_generator_)};
    const auto interactions{create_points(num_of_interactions)};
    const auto interaction_mark_position_pairs{convert_to_mark_position_pairs(interactions)};

    const auto connections{generate_connections(vertices, interactions)};
    const auto simplices{create_simplices_from_connections(connections)};
    FiniteNetwork network{max_dimension(), vertex_ids, simplices};

    return network;
}

PointList FiniteHypergraphModel::create_points(const size_t num_of_nodes) const
{
    const auto marks{generate_marks(num_of_nodes)};
    const auto positions{generate_positions(num_of_nodes)};

    std::uniform_real_distribution<> uniform_distribution_u(0., 1.);

    PointList nodes{};
    nodes.reserve(num_of_nodes);
    for (auto index{0U}; index < num_of_nodes; ++index)
    {
        nodes.push_back(Point(index, marks[index], positions[index]));
    }

    return nodes;
}

PositionList FiniteHypergraphModel::generate_positions(const size_t num_of_nodes) const
{
    std::uniform_real_distribution<Position> uniform_distribution(
        -TORUS_SIZE / 2., +TORUS_SIZE / 2.);

    PositionList positions(num_of_nodes, 0.0);
    for (auto i{0U}; i < num_of_nodes; ++i)
    {
        positions.at(i) = uniform_distribution(random_number_generator_);
    }

    return positions;
}

float FiniteHypergraphModel::distance(const Point &first, const Point &second) const
{
    const auto distance_inside{fabs(first.position() - second.position())};
    const auto distance{
        distance_inside < 0.5 * TORUS_SIZE ? distance_inside : TORUS_SIZE - distance_inside};

    return distance;
}