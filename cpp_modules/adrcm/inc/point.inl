#ifndef _POINT_INL_
#define _POINT_INL_

#include <cmath>

Point::Point(const double birth_time, const Position &position)
    : birth_time_(birth_time), position_(position)
{
}

inline double Point::distance(const Point &other, const double torus_size) const
{
    return dimension() == 1 ? distance_1d(other, torus_size) : distance_general(other, torus_size);
}

inline double Point::distance_general(const Point &other, const double torus_size) const
{
    const auto this_position{position()};
    const auto that_position{other.position()};

    auto sum_of_squared_coordinate_distances{0.};
    for (auto index{0U}; index < dimension(); ++index)
    {
        const auto current_dimension_distance_inside{this_position[index] - that_position[index]};
        const auto current_dimension_distance{
            current_dimension_distance_inside < 0.5 * torus_size ? current_dimension_distance_inside : torus_size - current_dimension_distance_inside};
        sum_of_squared_coordinate_distances += std::pow(
            current_dimension_distance,
            2U);
    }
    return std::pow(sum_of_squared_coordinate_distances, 0.5);
}

inline double Point::distance_1d(const Point &other, const double torus_size) const
{
    const auto this_position{position()};
    const auto that_position{other.position()};

    const auto distance_inside{fabs(this_position[0] - that_position[0])};

    const auto distance{
        distance_inside < 0.5 * torus_size ? distance_inside : torus_size - distance_inside};

    return distance;
}

inline const double &Point::birth_time() const
{
    return birth_time_;
}

inline const Position &Point::position() const
{
    return position_;
}

inline uint32_t Point::dimension() const
{
    return position_.size();
}

#endif