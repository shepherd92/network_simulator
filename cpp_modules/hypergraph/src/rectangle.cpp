#include <cassert>
#include <cmath>
#include <execution>
#include <iostream>
#include <mutex>

#include "point.h"
#include "rectangle.h"

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
      exponent_(other.exponent_),
      bottom_to_minus_exponent_(other.bottom_to_minus_exponent_),
      points_(std::move(other.points_))
{
}

void Rectangle::add_point(const Point &point)
{
    std::lock_guard<std::mutex> lock_guard(mutex);
    points_.push_back(point);
}

void Rectangle::set_exponent(const float exponent)
{
    std::lock_guard<std::mutex> lock_guard(mutex);
    exponent_ = exponent;
    bottom_to_minus_exponent_ = std::pow(bottom_mark_, exponent_);
}

bool Rectangle::contains(const Point &point) const
{
    return point.position() >= left() &&
           point.position() <= right() &&
           point.mark() >= bottom() &&
           point.mark() <= top();
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

const float &Rectangle::exponent() const
{
    return exponent_;
}

const float &Rectangle::bottom_to_minus_exponent() const
{
    return bottom_to_minus_exponent_;
}

const PointList &Rectangle::points() const
{
    return points_;
}

bool Rectangle::can_connect(
    const Rectangle &other,
    const float beta,
    const double torus_size) const
{
    // assumption: no rectangles can wrap around the torus boundaries
    if (left() <= other.right() && other.left() <= right())
    {
        // the distance is 0 if the intervals overlapif (distance < max_distance)
        return true;
    }

    const auto max_distance{
        beta * std::pow(bottom(), -exponent()) *
        std::pow(other.bottom(), -other.exponent())};

    const auto distance_inside{right() < other.left()
                                   ? other.left() - right()
                                   : left() - other.right()};
    const auto distance{
        distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside};

    return distance < max_distance;
}

ConnectionList Rectangle::calc_connected_point_pairs(
    const Rectangle &other,
    const float beta,
    const float torus_size) const
{
    ConnectionList connections{};
    if (!can_connect(other, beta, torus_size))
    {
        return connections;
    }

    std::for_each(
        std::execution::par_unseq,
        points().begin(),
        points().end(),
        [&](auto &&first_point)
        {
            const auto partial_result{beta * first_point.mark_to_gamma()};
            for (const auto &second_point : other.points())
            {
                const auto max_distance{partial_result * second_point.mark_to_gamma()};
                const auto distance{first_point.distance(second_point, torus_size)};
                if (distance < max_distance)
                {
                    std::lock_guard<std::mutex> lock_guard(mutex);
                    connections.push_back(std::pair(first_point.id(), second_point.id()));
                }
            }
        });

    return connections;
}
