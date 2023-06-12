#ifndef _ADRCM_MODEL_INL_
#define _ADRCM_MODEL_INL_

#include <cmath>

AdrcmModel::Vertex::Vertex(const vertex_id id, const float birth_time, const float position)
    : id_(id), birth_time_(birth_time), position_(position)
{
}

inline double AdrcmModel::Vertex::distance(const Vertex &other) const
{
    return std::abs(position() - other.position());
}

inline double AdrcmModel::Vertex::torus_distance(const Vertex &other, const double torus_size) const
{
    const auto distance_inside{distance(other)};
    assert(distance_inside < torus_size);
    return distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside;
}

inline const vertex_id &AdrcmModel::Vertex::id() const
{
    return id_;
}

inline const float &AdrcmModel::Vertex::birth_time() const
{
    return birth_time_;
}

inline const float &AdrcmModel::Vertex::position() const
{
    return position_;
}

inline bool AdrcmModel::Vertex::operator<(const Vertex &other) const
{
    return birth_time_ < other.birth_time();
}

inline double AdrcmModel::profile_function(const double argument) const
{
    return argument <= parameters.alpha ? 1. / (2. * parameters.alpha) : 0.;
}

#endif