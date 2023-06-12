#ifndef _ADRCM_MODEL_INL_
#define _ADRCM_MODEL_INL_

#include <cmath>

AdrcmModel::Point::Point(const vertex_id id, const float birth_time, const float position)
    : id_(id), birth_time_(birth_time), position_(position)
{
}

inline auto AdrcmModel::Point::distance(const Point &other) const
{
    return std::abs(position() - other.position());
}

inline auto AdrcmModel::Point::torus_distance(const Point &other, const double torus_size) const
{
    const auto distance_inside{distance(other)};
    assert(distance_inside < torus_size);
    return distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside;
}

inline auto AdrcmModel::Point::id() const
{
    return id_;
}

inline const auto &AdrcmModel::Point::birth_time() const
{
    return birth_time_;
}

inline const auto &AdrcmModel::Point::position() const
{
    return position_;
}

inline auto AdrcmModel::Point::operator<(const Point &other) const
{
    return birth_time_ < other.birth_time();
}

inline auto AdrcmModel::profile_function(const double argument) const
{
    return argument <= parameters.alpha ? 1. / (2. * parameters.alpha) : 0.;
}

#endif