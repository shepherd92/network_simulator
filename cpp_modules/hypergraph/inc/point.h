#ifndef _POINT_H_
#define _POINT_H_

#include "typedefs.h"

class Point
{
public:
    Point(const Id id, const Mark mark, const Position position, const float minus_exponent);
    bool connects(const Point &other, const float beta, const float torus_size) const;
    float distance(const Point &other, const float torus_size) const;

    const Id &id() const;
    const Mark &mark() const;
    const Position &position() const;
    const float &mark_to_gamma() const;

private:
    const Id id_;
    const Mark mark_;
    const Position position_;
    const float mark_to_gamma_;
};

#endif