#ifndef _RECTANGLE_H_
#define _RECTANGLE_H_

#include <cstdlib>
#include <vector>

#include "point.h"

class Rectangle
{
public:
    Rectangle(
        const float bottom_mark,
        const float top_mark,
        const float left_position,
        const float right_position);

    const float &top() const;
    const float &bottom() const;
    const float &left() const;
    const float &right() const;
    const float &exponent() const;
    const float &bottom_to_minus_exponent() const;

    void add_point(const Point &points);
    void set_exponent(const float exponent);

    const std::vector<Point> &points() const;

private:
    bool contains(const Point &point) const;

    const float bottom_mark_;
    const float top_mark_;
    const float left_position_;
    const float right_position_;
    float exponent_;
    float bottom_to_minus_exponent_;

    std::vector<Point> points_;
};

bool can_connect(
    const Rectangle &first,
    const Rectangle &second,
    const float beta,
    const double torus_size,
    const bool is_finite);

std::vector<std::pair<int, int>> calc_connected_point_pairs(
    const std::vector<Point> &points_1,
    const std::vector<Point> &points_2,
    const float b,
    const float torus_size,
    const bool is_finite);

#endif