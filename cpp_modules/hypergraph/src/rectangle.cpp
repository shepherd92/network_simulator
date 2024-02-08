#include <cassert>
#include <cmath>

#include "point.h"
#include "rectangle.h"

Rectangle::Rectangle(
    const float bottom_mark,
    const float top_mark,
    const float left_position,
    const float right_position)
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

const float &Rectangle::bottom() const
{
    return bottom_mark_;
}

const float &Rectangle::top() const
{
    return top_mark_;
}

const float &Rectangle::left() const
{
    return left_position_;
}

const float &Rectangle::right() const
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

bool can_connect(
    const Rectangle &first,
    const Rectangle &second,
    const float beta,
    const double torus_size,
    const bool is_finite)
{
    // assumption: no rectangles can wrap around the torus boundaries
    if (first.left() <= second.right() && second.left() <= first.right())
    {
        // the distance is 0 if the intervals overlap
        return true;
    }

    const auto max_distance{
        beta * std::pow(first.bottom(), -first.exponent()) *
        std::pow(second.bottom(), -second.exponent())};

    if (is_finite)
    {
        const auto distance_inside{first.right() < second.left()
                                       ? second.left() - first.right()
                                       : first.left() - second.right()};
        const auto distance{
            distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside};
        return distance < max_distance;
    }
    else
    {
        const auto distance{first.right() < second.left()
                                ? second.left() - first.right()
                                : first.left() - second.right()};
        return distance < max_distance;
    }
}

std::vector<std::pair<int, int>> calc_connected_point_pairs(
    const std::vector<Point> &points_1,
    const std::vector<Point> &points_2,
    const float b,
    const float torus_size,
    const bool is_finite)
{
    std::vector<std::pair<int, int>> connections{};
    for (const auto &first_point : points_1)
    {
        const auto partial_result{b * first_point.mark_to_gamma()};
        for (const auto &second_point : points_2)
        {
            const auto max_distance{partial_result * second_point.mark_to_gamma()};
            const auto distance{is_finite
                                    ? first_point.torus_distance(second_point, torus_size)
                                    : first_point.distance(second_point)};
            if (distance < max_distance)
            {
                connections.push_back(std::pair(first_point.id(), second_point.id()));
            }
        }
    }
    return connections;
}
