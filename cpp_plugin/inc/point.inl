#ifndef _POINT_INL_
#define _POINT_INL_

#include <cmath>

Point::Point(const Mark mark, const Position position, const PointId id)
    : mark_(mark), position_(position), id_(id)
{
}

const PointId &Point::id() const
{
    return id_;
}

const Mark &Point::mark() const
{
    return mark_;
}

const Position &Point::position() const
{
    return position_;
}

void Point::set_id(const PointId id)
{
    id_ = id;
}

void Point::set_mark(const Mark mark)
{
    mark_ = mark;
}

#endif