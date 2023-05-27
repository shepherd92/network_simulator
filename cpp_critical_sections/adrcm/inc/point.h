#ifndef _POINT_H_
#define _POINT_H_

#include <cstdlib>

class Point
{
public:
    Point(const double birth_time, const double position);

    inline double distance(const Point &other) const;
    inline double torus_distance(const Point &other, const double torus_size) const;

    inline const double &birth_time() const;
    inline const double &position() const;

private:
    const double birth_time_;
    const double position_;
};

#include "point.inl"

#endif