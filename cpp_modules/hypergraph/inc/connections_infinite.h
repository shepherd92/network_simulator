#ifndef _CONNECTIONS_INFINITE_H_
#define _CONNECTIONS_INFINITE_H_

#include <vector>

class Point;

std::vector<Point> create_interactions_infinite(
    const double u,
    const double beta,
    const double gamma,
    const double gamma_prime,
    const double interaction_intensity);

std::vector<Point> create_vertices_infinite(
    const std::vector<Point> &interactions,
    const double beta,
    const double gamma,
    const double gamma_prime);

double determine_infinite_network_size(
    const std::vector<Point> &interactions,
    const std::vector<Point> &vertices);

std::vector<float> generate_positions_infinite(
    const std::vector<float> &marks,
    const float beta_x_u_to_gamma,
    const float exponent);

#endif