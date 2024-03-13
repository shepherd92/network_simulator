#ifndef _RECTANGLE_H_
#define _RECTANGLE_H_

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

    ConnectionList calc_connected_point_pairs(const Rectangle &other, const float beta) const;

    const Mark &top() const;
    const Mark &bottom() const;
    const Position &left() const;
    const Position &right() const;
    const float &exponent() const;
    const float &bottom_to_minus_exponent() const;

    void add_point(const Point &points);
    void set_exponent(const float exponent);

    const PointList &points() const;

private:
    bool contains(const Point &point) const;

    const Mark bottom_mark_;
    const Mark top_mark_;
    const Position left_position_;
    const Position right_position_;
    float exponent_;
    float bottom_to_minus_exponent_;

    mutable std::mutex mutex;
    PointList points_;
};

void fill_rectangles(RectangleList &rectangles, const PointList &points);

#endif