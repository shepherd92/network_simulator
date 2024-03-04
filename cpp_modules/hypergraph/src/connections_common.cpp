#include <algorithm>
#include <execution>
#include <iostream>
#include <mutex>
#include <random>
#include <vector>

#include "connections_common.h"
#include "globals.h"
#include "point.h"
#include "rectangle.h"
#include "tools.h"

MarkList generate_marks(const size_t num_nodes, const Mark min_mark)
{
    std::uniform_real_distribution<Mark> uniform_distribution(0., 1.);
    MarkList marks(num_nodes, 0.0);

    for (auto i{0U}; i < num_nodes; ++i)
    {
        marks.at(i) = std::max(min_mark, uniform_distribution(random_number_generator));
    }

    std::sort(marks.begin(), marks.end());

    return marks;
}

RectangleList create_rectangles(const PointList &points, const double gamma)
{
    const auto n_points{points.size()};
    const auto torus_size{determine_network_size(points)};
    RectangleList rectangles{};
    const auto points_per_rectangle{std::max(100., std::pow(n_points, 0.5))};
    const auto area{points_per_rectangle / n_points * torus_size};

    auto bottom{0.0};
    while (bottom < 1.)
    {
        auto width{std::min(torus_size, std::pow(bottom, -gamma))};
        const auto height{std::min(area / width, 1. - bottom)};
        width = area / height; // adjust width to match the height if height is too large
        const auto rectangles_in_row{std::ceil(torus_size / width)};
        width = torus_size / rectangles_in_row; // adjust width to match the number of rectangles

        for (auto i{0U}; i < rectangles_in_row; ++i)
        {
            const auto left{i * width - torus_size / 2.};
            auto rectangle{Rectangle(bottom, bottom + height, left, left + width)};
            rectangle.set_exponent(gamma);
            rectangles.push_back(std::move(rectangle));
        }
        bottom += height;
    }
    fill_rectangles(rectangles, points);
    return rectangles;
}

void fill_rectangles(RectangleList &rectangles, const PointList &points)
{
    // assume points are sorted with respect to marks
    std::for_each(
        execution_policy,
        points.begin(),
        points.end(),
        [&](auto &&point)
        {
            auto min_rectangle_index{0U};
            auto current_bottom{0.};
            for (auto i{min_rectangle_index}; i < rectangles.size(); ++i)
            {
                if (!is_close(current_bottom, rectangles[i].bottom()))
                {
                    // we are now in the next row from the bottom
                    min_rectangle_index = i;
                    current_bottom = rectangles[i].bottom();
                }
                // assume that the rectangles are sorted, no need to check the bottom
                if (point.mark() <= rectangles[i].top() &&
                    point.position() >= rectangles[i].left() &&
                    point.position() <= rectangles[i].right())
                {
                    rectangles[i].add_point(point);
                    break;
                }
            }
        });
}

ConnectionList generate_connections(
    const RectangleList &interaction_rectangles,
    const RectangleList &vertex_rectangles,
    const double beta,
    const double torus_size)
{
    ConnectionList connections{};
    std::mutex mutex;

    std::atomic<uint32_t> counter{0U};
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
                    const auto new_connections{interaction_rectangle.calc_connected_point_pairs(
                        vertex_rectangle,
                        beta,
                        torus_size)};
                    std::lock_guard<std::mutex> lock_guard(mutex);
                    connections.insert(connections.end(), new_connections.begin(), new_connections.end());
                });
            if (++counter % 10000 == 0)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rGenerating connections: " << counter << " / " << interaction_rectangles.size();
            }
        });
    std::cout << "\rGenerating connections: " << counter << " / " << interaction_rectangles.size() << std::flush;

    return connections;
}

double determine_network_size(const PointList &points)
{
    Position max_position{0.};
    for (const auto &point : points)
    {
        max_position = std::max(max_position, Position(fabs(point.position())));
    }
    return 2. * max_position;
}