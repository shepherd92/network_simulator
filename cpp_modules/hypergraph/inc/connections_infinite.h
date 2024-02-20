#ifndef _CONNECTIONS_INFINITE_H_
#define _CONNECTIONS_INFINITE_H_

#include "typedefs.h"

PointList create_interactions_infinite(
    const Mark u,
    const double beta,
    const double gamma,
    const double gamma_prime,
    const double interaction_intensity);

PointList create_vertices_infinite(
    const PointList &interactions,
    const double beta,
    const double gamma);

PointList create_points_in_neighborhood(
    const Point &point,
    const double beta,
    const double other_exponent);

PositionList generate_positions_infinite(
    const MarkList &marks,
    const float beta_x_u_to_gamma,
    const float exponent);

#endif