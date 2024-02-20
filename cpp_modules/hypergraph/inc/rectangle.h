#ifndef _RECTANGLE_H_
#define _RECTANGLE_H_

#include "typedefs.h"

class Rectangle
{
public:
    Rectangle(
        const Mark bottom_mark,
        const Mark top_mark,
        const Position left_position,
        const Position right_position);

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

    PointList points_;
};

bool can_connect(
    const Rectangle &first,
    const Rectangle &second,
    const float beta,
    const double torus_size);

std::vector<std::pair<int, int>> calc_connected_point_pairs(
    const PointList &points_1,
    const PointList &points_2,
    const float b,
    const float torus_size);

#endif