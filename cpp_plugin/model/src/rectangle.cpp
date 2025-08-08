#include <cassert>
#include <cmath>
#include <iostream>
#include <mutex>

#include "point.hpp"
#include "rectangle.hpp"
#include "tools.hpp"
#include "typedefs.hpp"

Rectangle::Rectangle(
    const Mark bottom_mark,
    const Mark top_mark,
    const Position left_position,
    const Position right_position)
    : bottom_mark_(bottom_mark),
      top_mark_(top_mark),
      left_position_(left_position),
      right_position_(right_position)
{
    assert(bottom() < top());
    assert(left() < right());
}

Rectangle::Rectangle(Rectangle &&other)
    : bottom_mark_(other.bottom_mark_),
      top_mark_(other.top_mark_),
      left_position_(other.left_position_),
      right_position_(other.right_position_),
      points_(std::move(other.points_))
{
}

void Rectangle::add_point(const Point &point)
{
    std::lock_guard<std::mutex> lock_guard(mutex);
    points_.push_back(point);
}

bool Rectangle::contains(const Point &point) const
{
    return point.position() >= left() &&
           point.position() <= right() &&
           point.mark() >= bottom() &&
           point.mark() <= top();
}

void Rectangle::transform_points(const std::function<void(Point &)> &lambda)
{
    std::for_each(
        execution_policy,
        points_.begin(),
        points_.end(),
        [&](auto &point)
        {
            lambda(point);
        });
}

const Mark &Rectangle::bottom() const
{
    return bottom_mark_;
}

const Mark &Rectangle::top() const
{
    return top_mark_;
}

const Position &Rectangle::left() const
{
    return left_position_;
}

const Position &Rectangle::right() const
{
    return right_position_;
}

const PointList &Rectangle::points() const
{
    return points_;
}

void Rectangle::set_top(const Mark mark)
{
    top_mark_ = mark;
}

void Rectangle::set_bottom(const Mark mark)
{
    bottom_mark_ = mark;
}

void fill_rectangles(RectangleList &rectangles, const PointList &points)
{
    // assume points are sorted with respect to marks
    std::for_each(
        execution_policy,
        points.begin(),
        points.end(),
        [&](auto &&point)
        {
            auto min_rectangle_index{0U};
            auto current_bottom{0.};
            for (auto i{min_rectangle_index}; i < rectangles.size(); ++i)
            {
                if (!is_close(current_bottom, rectangles[i].bottom()))
                {
                    // we are now in the next row from the bottom
                    min_rectangle_index = i;
                    current_bottom = rectangles[i].bottom();
                }
                // assume that the rectangles are sorted, no need to check the bottom
                if (point.mark() <= rectangles[i].top() &&
                    point.position() >= rectangles[i].left() &&
                    point.position() <= rectangles[i].right())
                {
                    rectangles[i].add_point(point);
                    break;
                }
            }
        });
}

RectangleList remove_empty_rectangles(RectangleList &rectangles)
{
    RectangleList result{};
    std::for_each(
        std::execution::seq,
        rectangles.begin(),
        rectangles.end(),
        [&](auto &rectangle)
        {
            if (!rectangle.points().empty())
            {
                result.push_back(std::move(rectangle));
            }
        });
    return result;
}

void transform_points(RectangleList &rectangles, const std::function<void(Point &)> &lambda)
{
    std::for_each(
        std::execution::seq,
        rectangles.begin(),
        rectangles.end(),
        [&](auto &rectangle)
        {
            rectangle.transform_points(lambda);
        });
}

void transform_rectangles(RectangleList &rectangles, const std::function<void(Rectangle &)> &lambda)
{
    std::for_each(
        execution_policy,
        rectangles.begin(),
        rectangles.end(),
        [&](auto &rectangle)
        {
            lambda(rectangle);
        });
}

ConnectionList all_point_pairs(const Rectangle &left, const Rectangle &right)
{
    ConnectionList connections{};
    for (const auto &left_point : left.points())
    {
        for (const auto &right_point : right.points())
        {
            connections.emplace_back(left_point.id(), right_point.id());
        }
    }
    return connections;
}

std::ostream &operator<<(std::ostream &os, const Rectangle &rectangle)
{
    os << "Rectangle (b,t,l,r): " << rectangle.bottom() << " " << rectangle.top() << " " << rectangle.left() << " " << rectangle.right();
    return os;
}