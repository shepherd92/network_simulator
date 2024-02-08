#ifndef _POINT_H_
#define _POINT_H_

#include <cstdlib>

class Point
{
public:
    Point(const int id, const float birth_time, const float position, const float minus_exponent);

    bool connects(const Point &other, const float torus_size, const float beta, const bool is_finite) const;

    float distance(const Point &other) const;
    float torus_distance(const Point &other, const float torus_size) const;

    const int &id() const;
    const float &mark() const;
    const float &position() const;
    const float &mark_to_gamma() const;

private:
    const int id_;
    const float mark_;
    const float position_;
    const float mark_to_gamma_;
};

#endif