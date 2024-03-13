#ifndef _POINT_H_
#define _POINT_H_

#include <cstdlib>

#include "typedefs.h"

class Point
{
public:
    Point(const Mark birth_time, const Position &position);

    inline double distance(const Point &other, const double torus_size) const;

    inline uint32_t dimension() const;
    inline const double &birth_time() const;
    inline const Position &position() const;

private:
    const Mark mark_;
    const Position position_;
};

#include "adrcm_point.inl"

#endif