#ifndef _CONNECTIONS_COMMON_H_
#define _CONNECTIONS_COMMON_H_

#include "typedefs.h"

MarkList generate_marks(const size_t num_nodes);

RectangleList create_rectangles(
    const PointList &points,
    const double gamma);

void fill_rectangles(
    RectangleList &rectangles,
    const PointList &points);

ConnectionList generate_connections(
    const RectangleList &interaction_rectangles,
    const RectangleList &vertex_rectangles,
    const double beta,
    const double torus_size);

double determine_network_size(const PointList &points);

#endif