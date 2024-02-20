#ifndef _POINT_H_
#define _POINT_H_

#include <cstdlib>

#include "typedefs.h"

class Point
{
public:
    Point(const double birth_time, const Position &position);

    inline double distance(const Point &other, const double torus_size) const;

    inline uint32_t dimension() const;
    inline const double &birth_time() const;
    inline const Position &position() const;

private:
    inline double distance_1d(const Point &other, const double torus_size) const;
    inline double distance_general(const Point &other, const double torus_size) const;

    const double birth_time_;
    const Position position_;
};

#include "point.inl"

#endif