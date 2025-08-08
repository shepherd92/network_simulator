#ifndef _RECTANGLE_HPP_
#define _RECTANGLE_HPP_

#include <mutex>

#include "typedefs.hpp"

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
    void set_top(const Mark mark);
    void set_bottom(const Mark mark);

    void add_point(const Point &points);
    void transform_points(const std::function<void(Point &)> &lambda);

    const PointList &points() const;

private:
    bool contains(const Point &point) const;
    friend std::ostream &operator<<(std::ostream &os, const Rectangle &rectangle);

    Mark bottom_mark_;
    Mark top_mark_;
    const Position left_position_;
    const Position right_position_;

    mutable std::mutex mutex;
    PointList points_;
};

void fill_rectangles(RectangleList &rectangles, const PointList &points);
RectangleList remove_empty_rectangles(RectangleList &rectangles);
void transform_points(RectangleList &rectangles, const std::function<void(Point &)> &lambda);
void transform_rectangles(RectangleList &rectangles, const std::function<void(Rectangle &)> &lambda);
ConnectionList all_point_pairs(const Rectangle &left, const Rectangle &right);

#endif