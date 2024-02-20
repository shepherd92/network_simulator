#include <cmath>

#include "point.h"
#include <iostream>

Point::Point(const Id id, const Mark mark, const Position position, const float minus_exponent)
    : id_(id), mark_(mark), position_(position), mark_to_gamma_(std::pow(mark_, -minus_exponent))
{
}

bool Point::connects(const Point &other, const float beta, const float torus_size) const
{
    const auto d{distance(other, torus_size)};
    const auto max_distance_of_connection{mark_to_gamma() * other.mark_to_gamma() * beta};
    return d < max_distance_of_connection;
}

float Point::distance(const Point &other, const float torus_size) const
{
    const auto distance_inside{fabs(position() - other.position())};

    const auto distance{
        distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside};

    return distance;
}

const Id &Point::id() const
{
    return id_;
}

const Mark &Point::mark() const
{
    return mark_;
}

const Position &Point::position() const
{
    return position_;
}

const float &Point::mark_to_gamma() const
{
    return mark_to_gamma_;
}
