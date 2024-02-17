#ifndef _CONNECTIONS_COMMON_H_
#define _CONNECTIONS_COMMON_H_

class Point;
class Rectangle;
class ModelParameters;

std::vector<float> generate_marks(const size_t num_nodes);
std::vector<Rectangle> create_rectangles(
    const double torus_size,
    const uint32_t n_points,
    const double gamma);
void fill_rectangles(
    std::vector<Rectangle> &rectangles,
    const std::vector<Point> &points);
std::vector<std::pair<int, int>> generate_connections(
    const std::vector<Rectangle> &interaction_rectangles,
    const std::vector<Rectangle> &vertex_rectangles,
    const bool is_finite,
    const ModelParameters &model_parameters);

#endif