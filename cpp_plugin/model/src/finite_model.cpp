#include "finite_model.h"
#include "point.h"

FiniteModel::FiniteModel(const uint32_t seed) : Model{seed}
{
}

PointList FiniteModel::create_points(const size_t num_of_points) const
{
    const auto marks{generate_marks(num_of_points)};
    const auto positions{generate_positions(num_of_points)};

    PointList points{};
    points.reserve(num_of_points);
    for (auto index{0}; index < static_cast<int32_t>(num_of_points); ++index)
    {
        points.emplace_back(Point{marks[index], positions[index], index});
    }

    return points;
}

PositionList FiniteModel::generate_positions(const size_t num_of_points) const
{
    std::uniform_real_distribution<Position> uniform_distribution(
        -torus_size() / 2., +torus_size() / 2.);

    PositionList positions(num_of_points, 0.0);
    for (auto i{0U}; i < num_of_points; ++i)
    {
        positions.at(i) = uniform_distribution(random_number_generator_);
    }

    return positions;
}

float FiniteModel::distance(const Point &first, const Point &second) const
{
    const auto distance_inside{fabs(first.position() - second.position())};
    const auto distance{
        distance_inside < 0.5 * torus_size() ? distance_inside : torus_size() - distance_inside};
    return distance;
}
