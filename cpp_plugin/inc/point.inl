#ifndef _POINT_INL_
#define _POINT_INL_

Point::Point(const PointId id, const Mark mark, const Position position)
    : id_(id), mark_(mark), position_(position)
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

#endif