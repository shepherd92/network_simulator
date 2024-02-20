#include <cassert>
#include <cmath>
#include <iostream>

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

void Rectangle::add_point(const Point &point)
{
    points_.push_back(point);
}

void Rectangle::set_exponent(const float exponent)
{
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

bool can_connect(
    const Rectangle &first,
    const Rectangle &second,
    const float beta,
    const double torus_size)
{
    // assumption: no rectangles can wrap around the torus boundaries
    if (first.left() <= second.right() && second.left() <= first.right())
    {
        // the distance is 0 if the intervals overlapif (distance < max_distance)
        return true;
    }

    const auto max_distance{
        beta * std::pow(first.bottom(), -first.exponent()) *
        std::pow(second.bottom(), -second.exponent())};

    const auto distance_inside{first.right() < second.left()
                                   ? second.left() - first.right()
                                   : first.left() - second.right()};
    const auto distance{
        distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside};

    return distance < max_distance;
}

std::vector<std::pair<int, int>> calc_connected_point_pairs(
    const PointList &points_1,
    const PointList &points_2,
    const float b,
    const float torus_size)
{
    std::vector<std::pair<int, int>> connections{};
    for (const auto &first_point : points_1)
    {
        const auto partial_result{b * first_point.mark_to_gamma()};
        for (const auto &second_point : points_2)
        {
            const auto max_distance{partial_result * second_point.mark_to_gamma()};
            const auto distance{first_point.distance(second_point, torus_size)};
            if (distance < max_distance)
            {
                connections.push_back(std::pair(first_point.id(), second_point.id()));
            }
        }
    }
    return connections;
}