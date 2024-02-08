#include <iostream>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "connection_generator.h"
#include "model_parameters.h"
#include "numpy_cpp_conversion.h"
#include "point.h"
#include "rectangle.h"

namespace py = pybind11;

bool is_close(const double first, const double second);

std::vector<Point> create_nodes(
    const size_t num_of_nodes,
    const double torus_size,
    const float exponent,
    const uint32_t seed);

std::vector<float> generate_marks(const size_t num_nodes, const uint32_t seed);
std::vector<float> generate_positions(
    const size_t num_nodes,
    const double torus_size,
    const uint32_t seed);

std::vector<std::pair<int, int>> generate_connections(
    const std::vector<Point> &nodes,
    const std::vector<Point> &interactions,
    const bool is_finite,
    const ModelParameters &model_parameters,
    const uint32_t seed);

std::vector<Rectangle> create_rectangles(
    const double torus_size,
    const uint32_t n_horizontal,
    const uint32_t n_vertical);

void fill_rectangles(
    std::vector<Rectangle> &rectangles,
    const std::vector<Point> &points,
    const uint32_t n_cols);

std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>> generate_finite_network_connections_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed)
{
    const auto RECTANGLE_GRID_ROWS{10U};
    const auto RECTANGLE_GRID_COLS{10U};
    const auto model_parameters{ModelParameters(model_parameters_input)};

    const auto nodes{create_nodes(model_parameters.num_of_nodes, model_parameters.torus_size, model_parameters.gamma, seed)};
    auto vertex_rectangles{create_rectangles(model_parameters.torus_size, RECTANGLE_GRID_ROWS, RECTANGLE_GRID_COLS)};
    for (auto &rectangle : vertex_rectangles)
    {
        rectangle.set_exponent(model_parameters.gamma);
    }
    fill_rectangles(vertex_rectangles, nodes, RECTANGLE_GRID_COLS);

    const auto interactions{create_nodes(model_parameters.num_of_interactions, model_parameters.torus_size, model_parameters.gamma_prime, seed + 2U)};
    auto interaction_rectangles{create_rectangles(model_parameters.torus_size, RECTANGLE_GRID_ROWS, RECTANGLE_GRID_COLS)};
    for (auto &rectangle : interaction_rectangles)
    {
        rectangle.set_exponent(model_parameters.gamma_prime);
    }
    fill_rectangles(interaction_rectangles, interactions, RECTANGLE_GRID_COLS);

    const auto connections{generate_connections(nodes, interactions, true, model_parameters, seed)};

    const auto return_value_1{vector_of_pairs_to_numpy<int>(connections)};

    std::vector<std::pair<float, float>> node_mark_position_pairs{};
    node_mark_position_pairs.reserve(model_parameters.num_of_nodes);
    for (const auto &node : nodes)
    {
        node_mark_position_pairs.push_back(std::make_pair(node.mark(), node.position()));
    }
    const auto return_value_2{vector_of_pairs_to_numpy<float>(node_mark_position_pairs)};

    std::vector<std::pair<float, float>> interaction_mark_position_pairs{};
    interaction_mark_position_pairs.reserve(model_parameters.num_of_interactions);
    for (const auto &interaction : interactions)
    {
        interaction_mark_position_pairs.push_back(std::make_pair(interaction.mark(), interaction.position()));
    }
    const auto return_value_3{vector_of_pairs_to_numpy<float>(interaction_mark_position_pairs)};

    return std::make_tuple(std::move(return_value_1), std::move(return_value_2), std::move(return_value_3));
}

std::vector<Point> create_nodes(const size_t num_of_nodes, const double torus_size, const float exponent, const uint32_t seed)
{
    std::vector<uint32_t> node_ids(num_of_nodes, 0U);
    const auto marks{generate_marks(num_of_nodes, seed)};
    const auto positions{generate_positions(num_of_nodes, torus_size, seed + 1U)};

    std::uniform_real_distribution<> uniform_distribution_u(0., 1.);
    std::mt19937 random_number_generator{seed};

    std::vector<Point> nodes{};
    nodes.reserve(num_of_nodes);
    for (auto index{0U}; index < num_of_nodes; ++index)
    {
        nodes.push_back(Point(index, marks[index], positions[index], exponent));
    }

    return nodes;
}

std::vector<float> generate_marks(const size_t num_nodes, const uint32_t seed)
{
    std::mt19937 random_number_generator{seed};
    std::uniform_real_distribution<float> uniform_distribution(0., 1.);
    std::vector<float> marks(num_nodes, 0.0);

    for (auto i{0U}; i < num_nodes; ++i)
    {
        marks.at(i) = uniform_distribution(random_number_generator);
    }

    std::sort(marks.begin(), marks.end());

    return marks;
}

std::vector<float> generate_positions(
    const size_t num_of_nodes,
    const double torus_size,
    const uint32_t seed)
{
    std::mt19937 random_number_generator{seed};
    std::uniform_real_distribution<float> uniform_distribution(
        -torus_size / 2., +torus_size / 2.);

    std::vector<float> positions(num_of_nodes, 0.0);
    for (auto i{0U}; i < num_of_nodes; ++i)
    {
        positions.at(i) = uniform_distribution(random_number_generator);
    }

    return positions;
}

std::vector<std::pair<int, int>> generate_connections(
    const std::vector<Rectangle> &vertex_rectangles,
    const std::vector<Rectangle> &interaction_rectangles,
    const bool is_finite,
    const ModelParameters &model_parameters,
    const uint32_t seed)
{
    static std::mt19937 random_number_generator{seed};
    std::vector<std::pair<int, int>> connections{};

    auto counter{0U};
    for (const auto &interaction_rectangle : interaction_rectangles)
    {
        for (const auto &vertex_rectangle : vertex_rectangles)
        {
            if (can_connect(
                    interaction_rectangle,
                    vertex_rectangle,
                    model_parameters.beta,
                    model_parameters.torus_size,
                    is_finite))
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
        std::cout << "\rGenerating connections: " << counter << " / " << interaction_rectangles.size() << std::flush;
        ++counter;
    }
    return connections;
}

std::vector<Rectangle> create_rectangles(const double torus_size, const uint32_t n_cols, const uint32_t n_rows)
{
    const auto column_width{torus_size / n_cols};
    const auto row_height{1. / n_rows};

    std::vector<Rectangle> rectangles{};
    rectangles.reserve(n_rows * n_cols);
    for (auto row{0U}; row < n_rows; ++row)
    {
        for (auto column{0U}; column < n_cols; ++column)
        {
            rectangles.push_back(Rectangle(
                row * row_height,
                (row + 1U) * row_height,
                column * column_width - torus_size / 2,
                (column + 1U) * column_width - torus_size / 2));
        }
    }
    return rectangles;
}

void fill_rectangles(
    std::vector<Rectangle> &rectangles,
    const std::vector<Point> &points,
    const uint32_t n_cols)
{
    const auto row_height{rectangles[0].top() - rectangles[0].bottom()};
    const auto column_width{rectangles[0].right() - rectangles[0].left()};
    const auto min_position{rectangles[0].left()};

    for (const auto &point : points)
    {
        const auto row_index{static_cast<uint32_t>(point.mark() / row_height)};
        const auto column_index{static_cast<uint32_t>((point.position() - min_position) / column_width)};
        rectangles[row_index * n_cols + column_index].add_point(point);
    }
}
