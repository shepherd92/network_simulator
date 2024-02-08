#include <cmath>

#include "point.h"

Point::Point(const int id, const float mark, const float position, const float minus_exponent)
    : id_(id), mark_(mark), position_(position), mark_to_gamma_(std::pow(mark_, -minus_exponent))
{
}

bool Point::connects(const Point &other, const float torus_size, const float beta, const bool is_finite) const
{
    const auto d{is_finite ? torus_distance(other, torus_size) : distance(other)};
    const auto max_distance_of_connection{mark_to_gamma() * other.mark_to_gamma() * beta};
    return d < max_distance_of_connection;
}

float Point::distance(const Point &other) const
{
    return abs(position() - other.position());
}

float Point::torus_distance(const Point &other, const float torus_size) const
{
    const auto distance_inside{abs(position() - other.position())};

    const auto distance{
        distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside};

    return distance;
}

const int &Point::id() const
{
    return id_;
}

const float &Point::mark() const
{
    return mark_;
}

const float &Point::position() const
{
    return position_;
}

const float &Point::mark_to_gamma() const
{
    return mark_to_gamma_;
}
