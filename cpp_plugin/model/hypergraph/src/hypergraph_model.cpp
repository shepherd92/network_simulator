#include <cassert>
#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <vector>

#include "hypergraph_model.h"
#include "point.h"
#include "rectangle.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"
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

HypergraphModel::Parameters HypergraphModel::parameters() const
{
    return parameters_;
}

ConnectionList HypergraphModel::generate_connections(
    const PointList &vertices,
    const PointList &interactions) const
{
    const auto vertex_rectangles{create_transformed_filled_rectangles(vertices, -gamma())};
    const auto interaction_rectangles{create_transformed_filled_rectangles(interactions, -gamma_prime())};

    std::mutex mutex{};
    ConnectionList all_connections{};
    std::atomic_size_t progress{0};
    std::for_each(
        execution_policy,
        interaction_rectangles.begin(),
        interaction_rectangles.end(),
        [&](auto &&interaction_rectangle)
        {
            ConnectionList connections{};
            std::for_each(
                std::execution::seq,
                vertex_rectangles.begin(),
                vertex_rectangles.end(),
                [&](auto &&vertex_rectangle)
                {
                    const auto new_connections{calc_connected_point_pairs(vertex_rectangle, interaction_rectangle)};
                    connections.insert(connections.end(), new_connections.begin(), new_connections.end());
                });

            std::lock_guard<std::mutex> lock{mutex};
            all_connections.insert(all_connections.end(), connections.begin(), connections.end());
            log_progress(++progress, interaction_rectangles.size(), 10, "Generating connections");
        });
    return all_connections;
}

RectangleList HypergraphModel::create_transformed_filled_rectangles(
    const PointList &points,
    const float exponent) const
{
    const auto sqrt_beta{std::pow(beta(), 0.5)};
    auto rectangles{create_rectangles(points, exponent)};
    fill_rectangles(rectangles, points);
    rectangles = remove_empty_rectangles(rectangles);
    transform_points(rectangles,
                     [&](Point &point)
                     { point.set_mark(sqrt_beta * std::pow(point.mark(), exponent)); });
    transform_rectangles(rectangles, [&](Rectangle &rectangle)
                         { rectangle.set_top(sqrt_beta * std::pow(rectangle.top(), exponent));
                           rectangle.set_bottom(sqrt_beta * std::pow(rectangle.bottom(), exponent)); });
    return rectangles;
}

RectangleList HypergraphModel::create_rectangles(const PointList &points_in, const float exponent) const
{
    const auto n_points{points_in.size()};
    const auto space_size{std::max(determine_space_size(points_in), 1e-9F)};
    RectangleList rectangles{};
    if (n_points < 1000U)
    {
        rectangles.emplace_back(Rectangle{0.F, 1.F, -space_size / 2.F, space_size / 2.F});
        return rectangles;
    }
    const auto points_per_rectangle{std::max(100.F, std::pow(static_cast<float>(n_points), 0.5F))};
    const auto area{points_per_rectangle / n_points * space_size};
    const auto sqrt_beta{std::pow(beta(), 0.5F)};

    auto bottom{0.0F};
    while (bottom < 1.F)
    {
        auto width{std::min(space_size, sqrt_beta * std::pow(bottom, exponent))};
        const auto height{std::max(std::min(area / width, 1.F - bottom), 1e-6F)};
        const auto rectangles_in_row{std::ceil(1.F / width)};
        width = space_size / rectangles_in_row; // adjust width to match the number of rectangles

        for (auto i{0U}; i < rectangles_in_row; ++i)
        {
            const auto left{i * width - space_size / 2.F};
            rectangles.emplace_back(Rectangle(bottom, bottom + height, left, left + width));
        }
        bottom += height;
    }
    return rectangles;
}

ConnectionList HypergraphModel::calc_connected_point_pairs(
    const Rectangle &vertex_rectangle,
    const Rectangle &interaction_rectangle) const
{
    if (!rectangle_points_may_connect(vertex_rectangle, interaction_rectangle))
    {
        return ConnectionList{};
    }

    if (rectangle_points_surely_connect(vertex_rectangle, interaction_rectangle))
    {
        return all_point_pairs(vertex_rectangle, interaction_rectangle);
    }

    ConnectionList connections{};
    std::for_each(
        std::execution::seq,
        vertex_rectangle.points().begin(),
        vertex_rectangle.points().end(),
        [&](auto &&vertex)
        {
            std::for_each(
                std::execution::seq,
                interaction_rectangle.points().begin(),
                interaction_rectangle.points().end(),
                [&](auto &&interaction)
                {
                    if (connects(vertex, interaction))
                    {
                        connections.emplace_back(vertex.id(), interaction.id());
                    }
                });
        });

    return connections;
}

bool HypergraphModel::rectangle_points_may_connect(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const
{
    // assumption: no rectangles wrap around the torus boundaries
    if (vertex_rectangle.left() <= interaction_rectangle.right() && interaction_rectangle.left() <= vertex_rectangle.right())
    {
        // the distance is 0 if the intervals overlap
        return true;
    }

    return connects(Point{vertex_rectangle.bottom(), vertex_rectangle.right()},
                    Point{interaction_rectangle.bottom(), interaction_rectangle.left()}) ||
           connects(Point{vertex_rectangle.bottom(), vertex_rectangle.right()},
                    Point{interaction_rectangle.bottom(), interaction_rectangle.left()});
}

SimplexList HypergraphModel::create_interaction_simplices_from_connections(
    const PointIdList &interactions,
    const ConnectionList &connections) const
{
    // initialize the map with empty vectors
    std::map<PointId, std::vector<PointId>> interactions_map{};
    for (const auto &interaction : interactions)
    {
        interactions_map[interaction];
    }

    // fill the map with the connected vertices
    for (const auto &pair : connections)
    {
        interactions_map[pair.second].push_back(pair.first);
    }

    // create simplices from the map
    std::vector<Simplex> simplices{};
    for (const auto &interaction : interactions_map)
    {
        simplices.push_back(Simplex{interaction.second});
    }

    return SimplexList{simplices};
}
