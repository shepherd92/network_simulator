#ifndef _ADRCM_MODEL_INL_
#define _ADRCM_MODEL_INL_

#include <cmath>

inline float AdrcmModel::Point::distance(const Point &other) const
{
    return std::abs(position() - other.position());
}

inline float AdrcmModel::Point::torus_distance(const Point &other, const double torus_size) const
{
    const auto distance_inside{distance(other)};
    assert(distance_inside < torus_size);
    return distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside;
}

inline vertex_id AdrcmModel::Point::id() const
{
    return id_;
}

inline const float &AdrcmModel::Point::birth_time() const
{
    return birth_time_;
}

inline const float &AdrcmModel::Point::position() const
{
    return position_;
}

inline bool AdrcmModel::Point::operator<(const Point &other) const
{
    return birth_time_ < other.birth_time();
}

inline double AdrcmModel::profile_function(const double argument) const
{
    return argument <= parameters.alpha ? 1. / (2. * parameters.alpha) : 0.;
}

#endif