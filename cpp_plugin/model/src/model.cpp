#include "model.hpp"
#include "point.hpp"

Model::Model(const uint32_t seed)
    : random_number_generator_{seed}
{
}

MarkList Model::generate_marks(const size_t num_nodes, const Mark min_mark, const Mark max_mark) const
{
    std::uniform_real_distribution<Mark> uniform_distribution(0., max_mark);
    MarkList marks(num_nodes, 0.0);

    for (auto i{0U}; i < num_nodes; ++i)
    {
        marks.at(i) = std::max(min_mark, uniform_distribution(random_number_generator_));
    }

    std::sort(marks.begin(), marks.end());

    return marks;
}

float Model::determine_space_size(const PointList &points) const
{
    Position max_position{0.};
    for (const auto &point : points)
    {
        max_position = std::max(max_position, Position(fabs(point.position())));
    }
    return 2.F * max_position;
}
