#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "connections_common.h"
#include "globals.h"
#include "model_parameters.h"
#include "rectangle.h"
#include "tools.h"

std::vector<float> generate_marks(const size_t num_nodes)
{
    std::uniform_real_distribution<float> uniform_distribution(0., 1.);
    std::vector<float> marks(num_nodes, 0.0);

    for (auto i{0U}; i < num_nodes; ++i)
    {
        marks.at(i) = uniform_distribution(random_number_generator);
    }

    std::sort(marks.begin(), marks.end());

    return marks;
}

std::vector<Rectangle> create_rectangles(
    const double torus_size,
    const uint32_t n_points,
    const double gamma)
{
    std::vector<Rectangle> rectangles{};
    const auto points_per_rectangle{std::max(100., std::pow(n_points, 0.5))};
    const auto area{points_per_rectangle * torus_size / n_points};

    auto bottom{0.0};
    while (bottom < 1.)
    {
        const auto ideal_width{std::min(torus_size, std::pow(bottom, -gamma))};
        const auto rectangles_in_row{std::ceil(torus_size / ideal_width)};
        const auto width{torus_size / rectangles_in_row};
        const auto height{area / width};
        const auto top{std::min(1., bottom + height)};
        for (auto i{0U}; i < rectangles_in_row; ++i)
        {
            const auto left{i * width - torus_size / 2.};
            const auto right{(i + 1U) * width - torus_size / 2.};
            // std::cout << bottom << " " << top << " " << left << " " << right << std::endl;
            auto rectangle{Rectangle(bottom, top, left, right)};
            rectangle.set_exponent(gamma);
            rectangles.push_back(rectangle);
        }
        bottom += height;
    }
    return rectangles;
}

void fill_rectangles(
    std::vector<Rectangle> &rectangles,
    const std::vector<Point> &points)
{
    // assume points are sorted with respect to marks
    for (const auto &point : points)
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
    }
}

std::vector<std::pair<int, int>> generate_connections(
    const std::vector<Rectangle> &interaction_rectangles,
    const std::vector<Rectangle> &vertex_rectangles,
    const bool is_finite,
    const ModelParameters &model_parameters)
{
    std::vector<std::pair<int, int>> connections{};

    auto counter{0U};
    for (const auto &interaction_rectangle : interaction_rectangles)
    {
        for (const auto &vertex_rectangle : vertex_rectangles)
        {
            const auto connections_possible{can_connect(
                interaction_rectangle,
                vertex_rectangle,
                model_parameters.beta,
                model_parameters.torus_size,
                is_finite)};

            if (connections_possible)
            {
                const auto new_connections{calc_connected_point_pairs(
                    interaction_rectangle.points(),
                    vertex_rectangle.points(),
                    model_parameters.beta,
                    model_parameters.torus_size,
                    is_finite)};
                connections.insert(connections.end(), new_connections.begin(), new_connections.end());
            }
        }
        ++counter;
        std::cout << "\rGenerating connections: " << counter << " / " << interaction_rectangles.size() << std::flush;
    }
    return connections;
}