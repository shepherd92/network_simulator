
#ifndef _POINT_INL_
#define _POINT_INL_

#include <cmath>

Point::Point(const double birth_time, const double position)
    : birth_time_(birth_time), position_(position)
{
}

inline double Point::distance(const Point &other) const
{
    return std::abs(position() - other.position());
}

inline double Point::torus_distance(const Point &other, const double torus_size) const
{
    const auto distance_inside{distance(other)};
    assert(distance_inside < torus_size);
    return distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside;
}

inline const double &Point::birth_time() const
{
    return birth_time_;
}

inline const double &Point::position() const
{
    return position_;
}

#endif