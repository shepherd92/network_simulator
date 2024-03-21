#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <vector>

#include "hypergraph_model.h"
#include "point.h"
#include "rectangle.h"
#include "simplex.h"
#include "typedefs.h"

HypergraphModel::Parameters::Parameters(const std::vector<double> &parameters)
{
    max_dimension = static_cast<Dimension>(parameters[0]);
    network_size = parameters[1];
    interaction_intensity = parameters[2];
    beta = parameters[3];
    gamma = parameters[4];
    gamma_prime = parameters[5];
}

HypergraphModel::HypergraphModel(const std::vector<double> &parameters)
    : parameters_{parameters}
{
}

ConnectionList HypergraphModel::generate_connections(
    const PointList &vertices,
    const PointList &interactions) const
{
    auto vertex_rectangles{create_rectangles(vertices, -gamma())};
    fill_rectangles(vertex_rectangles, vertices);
    transform_points(vertex_rectangles,
                     [&](Point &point)
                     { point.set_mark(beta() * std::pow(point.mark(), -gamma())); });
    auto interaction_rectangles{create_rectangles(interactions, -gamma_prime())};
    fill_rectangles(interaction_rectangles, interactions);
    transform_points(interaction_rectangles,
                     [&](Point &point)
                     { point.set_mark(std::pow(point.mark(), -gamma_prime())); });

    ConnectionList connections{};
    std::for_each(
        std::execution::seq,
        interaction_rectangles.begin(),
        interaction_rectangles.end(),
        [&](auto &&interaction_rectangle)
        {
            std::for_each(
                std::execution::seq,
                vertex_rectangles.begin(),
                vertex_rectangles.end(),
                [&](auto &&vertex_rectangle)
                {
                    const auto new_connections{calc_connected_point_pairs(vertex_rectangle, interaction_rectangle)};
                    connections.insert(connections.end(), new_connections.begin(), new_connections.end());
                });
        });
    return connections;
}

RectangleList HypergraphModel::create_rectangles(const PointList &points, const float exponent) const
{
    const auto n_points{points.size()};
    RectangleList rectangles{};
    const auto points_per_rectangle{std::pow(n_points, 0.5)};
    const auto space_size{determine_space_size(points)};
    const auto area{points_per_rectangle / n_points * space_size};

    auto bottom{0.0};
    while (bottom < 1.)
    {
        auto width{std::min(1., std::pow(bottom, exponent))};
        const auto height{std::min(area / width, 1. - bottom)};
        width = area / height; // adjust width to match the height if height is too large
        const auto rectangles_in_row{std::ceil(1. / width)};
        width = space_size / rectangles_in_row; // adjust width to match the number of rectangles

        for (auto i{0U}; i < rectangles_in_row; ++i)
        {
            const auto left{i * width - space_size / 2.};
            auto rectangle{Rectangle(bottom, bottom + height, left, left + width)};
            rectangles.push_back(std::move(rectangle));
        }
        bottom += height;
    }
    return rectangles;
}

ConnectionList HypergraphModel::calc_connected_point_pairs(
    const Rectangle &vertex_rectangle,
    const Rectangle &interaction_rectangle) const
{
    if (!connects(vertex_rectangle, interaction_rectangle))
    {
        return ConnectionList{};
    }

    std::mutex mutex{};
    ConnectionList connections{};
    std::for_each(
        execution_policy,
        vertex_rectangle.points().begin(),
        vertex_rectangle.points().end(),
        [&](auto &&vertex)
        {
            for (const auto &interaction : interaction_rectangle.points())
            {
                if (connects(vertex, interaction))
                {
                    std::lock_guard<std::mutex> lock_guard(mutex);
                    connections.emplace_back(std::pair{vertex.id(), interaction.id()});
                }
            }
        });

    return connections;
}

bool HypergraphModel::connects(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const
{
    // assumption: no rectangles can wrap around the torus boundaries
    if (vertex_rectangle.left() <= interaction_rectangle.right() && interaction_rectangle.left() <= vertex_rectangle.right())
    {
        // the distance is 0 if the intervals overlap
        return true;
    }

    // add a dummy ID of 0 to the corner positions
    return connects(Point(0, beta() * std::pow(vertex_rectangle.bottom(), -gamma()), vertex_rectangle.right()),
                    Point(0, std::pow(interaction_rectangle.bottom(), -gamma_prime()), interaction_rectangle.left())) ||
           connects(Point(0, beta() * std::pow(vertex_rectangle.bottom(), -gamma()), vertex_rectangle.right()),
                    Point(0, std::pow(interaction_rectangle.bottom(), -gamma_prime()), interaction_rectangle.left()));
}

bool HypergraphModel::connects(const Point &vertex, const Point &interaction) const
{
    const auto dist{distance(vertex, interaction)};
    const auto max_distance_of_connection{vertex.mark() * interaction.mark()};
    return dist < max_distance_of_connection;
}

SimplexList HypergraphModel::create_simplices_from_connections(const ConnectionList &connections) const
{
    std::map<PointId, std::vector<PointId>> interactions;
    for (const auto &pair : connections)
    {
        interactions[pair.second].push_back(pair.first);
    }

    SimplexList simplices{};
    for (const auto &interaction : interactions)
    {
        simplices.push_back(Simplex{interaction.second});
    }

    return simplices;
}
