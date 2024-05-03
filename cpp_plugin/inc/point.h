#ifndef _POINT_H_
#define _POINT_H_

#include "typedefs.h"

class Point
{
public:
    inline Point(const Mark mark, const Position position, const PointId id = -1);

    inline void raise_mark_to_power(const float exponent);

    inline const PointId &id() const;
    inline const Mark &mark() const;
    inline const Position &position() const;

    inline void set_id(const PointId id);
    inline void set_mark(const Mark mark);

private:
    Mark mark_;
    Position position_;
    PointId id_;
};

#include "point.inl"

PointIdList convert_to_id_list(const PointList &points);
MarkList convert_to_mark_list(const PointList &points);
MarkPositionList convert_to_mark_position_pairs(const PointList &points);
void transform_points(PointList &points, const std::function<void(Point &)> &lambda);

#endif