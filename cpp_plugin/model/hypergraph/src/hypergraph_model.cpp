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

ConnectionList HypergraphModel::generate_connections(
    const PointList &vertices,
    const PointList &interactions) const
{
    const auto vertex_rectangles{create_transformed_filled_rectangles(vertices, -gamma())};
    const auto interaction_rectangles{create_transformed_filled_rectangles(interactions, -gamma_prime())};

    std::mutex mutex{};
    ConnectionList connections{};
    std::atomic_size_t progress{0};
    std::for_each(
        std::execution::seq,
        interaction_rectangles.begin(),
        interaction_rectangles.end(),
        [&](auto &&interaction_rectangle)
        {
            std::for_each(
                execution_policy,
                vertex_rectangles.begin(),
                vertex_rectangles.end(),
                [&](auto &&vertex_rectangle)
                {
                    const auto new_connections{calc_connected_point_pairs(vertex_rectangle, interaction_rectangle)};
                    std::lock_guard<std::mutex> lock{mutex};
                    connections.insert(connections.end(), new_connections.begin(), new_connections.end());
                });
            log_progress(++progress, interaction_rectangles.size(), 10, "Generating connections");
        });
    return connections;
}

RectangleList HypergraphModel::create_transformed_filled_rectangles(
    const PointList &points,
    const float exponent) const
{
    auto rectangles{create_rectangles(points, exponent)};
    fill_rectangles(rectangles, points);
    transform_points(rectangles,
                     [&](Point &point)
                     { point.set_mark(std::pow(beta(), 0.5) * std::pow(point.mark(), exponent)); });
    // transform_rectangles(rectangles, [&](Rectangle &rectangle)
    //                      { rectangle.set_top(std::pow(beta(), 0.5) * std::pow(rectangle.top(), exponent));
    //                        rectangle.set_bottom(std::pow(beta(), 0.5) * std::pow(rectangle.bottom(), exponent)); });
    return rectangles;
}

RectangleList HypergraphModel::create_rectangles(
    const PointList &points_in,
    const float exponent) const
{
    const auto n_points{points_in.size()};
    RectangleList rectangles{};
    // const auto points_per_rectangle{std::max(100., std::pow(n_points, 0.5))};
    const auto points_per_rectangle{n_points / 10.};
    const auto space_size{determine_space_size(points_in)};
    const auto area{points_per_rectangle / n_points * space_size};

    auto bottom{0.0};
    while (bottom < 1.)
    {
        auto width{std::min(1., std::pow(beta(), 0.5) * std::pow(bottom, exponent))};
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
    if (!rectangle_points_may_connect(vertex_rectangle, interaction_rectangle))
    {
        return ConnectionList{};
    }

    // if (rectangle_points_surely_connect(vertex_rectangle, interaction_rectangle))
    // {
    //     return all_point_pairs(vertex_rectangle, interaction_rectangle);
    // }

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

bool HypergraphModel::rectangle_points_surely_connect(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const
{
    const auto sqrt_beta{std::pow(beta(), 0.5F)};

    std::cout << interaction_rectangle << std::endl;

    const Point ver_top_left{sqrt_beta * std::pow(vertex_rectangle.top(), -gamma()), vertex_rectangle.left()};
    const Point ver_top_right{sqrt_beta * std::pow(vertex_rectangle.top(), -gamma()), vertex_rectangle.right()};
    const Point int_top_left{sqrt_beta * std::pow(interaction_rectangle.top(), -gamma_prime()), interaction_rectangle.left()};
    const Point int_top_right{sqrt_beta * std::pow(interaction_rectangle.top(), -gamma_prime()), interaction_rectangle.right()};

    const auto surely_connects{connects(ver_top_left, int_top_left) &&
                               connects(ver_top_left, int_top_right) &&
                               connects(ver_top_right, int_top_left) &&
                               connects(ver_top_right, int_top_right)};

    //////////////////////////////////////////////////////////////////////////
    if (surely_connects)
    {
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
                        if (!connects(vertex, interaction))
                        {
                            assert(false);
                        }
                    });
            });
    }

    //////////////////////////////////////////////////////////////////////////
    return surely_connects;
}

bool HypergraphModel::rectangle_points_may_connect(const Rectangle &vertex_rectangle, const Rectangle &interaction_rectangle) const
{
    // assumption: no rectangles wrap around the torus boundaries
    if (vertex_rectangle.left() <= interaction_rectangle.right() && interaction_rectangle.left() <= vertex_rectangle.right())
    {
        // the distance is 0 if the intervals overlap
        return true;
    }

    const auto sqrt_beta{std::pow(beta(), 0.5F)};

    return connects(Point{sqrt_beta * std::pow(vertex_rectangle.bottom(), -gamma()), vertex_rectangle.right()},
                    Point{sqrt_beta * std::pow(interaction_rectangle.bottom(), -gamma_prime()), interaction_rectangle.left()}) ||
           connects(Point{sqrt_beta * std::pow(vertex_rectangle.bottom(), -gamma()), vertex_rectangle.right()},
                    Point{sqrt_beta * std::pow(interaction_rectangle.bottom(), -gamma_prime()), interaction_rectangle.left()});
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