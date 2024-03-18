#ifndef _RECTANGLE_H_
#define _RECTANGLE_H_

#include <mutex>

#include "typedefs.h"

class Point;

class Rectangle
{
public:
    Rectangle(
        const Mark bottom_mark,
        const Mark top_mark,
        const Position left_position,
        const Position right_position);

    Rectangle(Rectangle &&other);

    const Mark &top() const;
    const Mark &bottom() const;
    const Position &left() const;
    const Position &right() const;

    void add_point(const Point &points);

    const PointList &points() const;

private:
    bool contains(const Point &point) const;

    const Mark bottom_mark_;
    const Mark top_mark_;
    const Position left_position_;
    const Position right_position_;

    mutable std::mutex mutex;
    PointList points_;
};

void fill_rectangles(RectangleList &rectangles, const PointList &points);

#endif