#include <cmath>

#include "hypergraph_model.h"
#include "neighborhood.h"

NeighborhoodPart::NeighborhoodPart(const float left, const float right)
    : left{left}, right{right}
{
}

Position NeighborhoodPart::left() const
{
    return left;
}

Position NeighborhoodPart::right() const
{
    return right;
}

void NeighborhoodPart::set_left(const Position value)
{
    left = value;
}

void NeighborhoodPart::set_right(const Position value)
{
    right = value;
}

Hyperbola::Hyperbola(
    const Position left,
    const Position right,
    const Position position,
    const Mark transformed_mark)
    : NeighborhoodPart{left, right}, position{position}, transformed_mark{transformed_mark}
{
}

float Hyperbola::operator()(const Position x, const float gamma) const
{
    if (x < left() || x > right())
    {
        return NAN;
    }
    if (!std::isfinite(x))
    {
        return 0.;
    }
    return std::pow(std::abs(position - x) / transformed_mark, -1.F / gamma);
}

PointList Hyperbola::create_points(const HypergraphModel::Parameters &parameters, std::mt19937 &rng) const
{
    // aliases
    const auto g{parameters.gamma};
    const auto uniform_exponent{g / (g - 1.F)};

    PointList vertices{};
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(parameters.network_size * integral(g))(rng)};
    vertices.reserve(num_of_vertices);

    std::uniform_real_distribution<float> uniform_01{};

    const auto uniform_lower_limit{std::pow(std::abs(right() - position), 1.F - 1.F / g)};
    const auto uniform_upper_limit{std::pow(std::abs(left() - position), 1.F - 1.F / g)};
    std::uniform_real_distribution<float> uniform{uniform_lower_limit, uniform_upper_limit};

    const auto transformation_coefficient{left() < position ? -1.F : 1.F};
    for (auto i{0U}; i < num_of_vertices; ++i)
    {
        // transform Z to be the position later
        const auto Z{uniform(rng)};
        const auto vertex_position{position + transformation_coefficient * std::pow(Z, uniform_exponent)};
        const auto max_mark{(*this)(vertex_position, g)};
        const auto mark{max_mark * uniform_01(rng)};
        vertices.emplace_back(Point{mark, vertex_position});
    }
    return vertices;
}

float Hyperbola::integral(const float gamma) const
{
    return gamma / (1. - gamma) *
           std::pow(transformed_mark, 1. / gamma) *
           std::abs((std::pow(std::abs(position - right()), -(1. - gamma) / gamma) -
                     std::pow(std::abs(position - left()), -(1. - gamma) / gamma)));
}

Center::Center(const float left, const float right)
    : NeighborhoodPart{left, right}
{
}

PointList Center::create_points(const HypergraphModel::Parameters &parameters, std::mt19937 &rng) const
{
    PointList vertices{};
    const auto num_of_vertices{std::poisson_distribution<uint32_t>(parameters.network_size * (right() - left()))(rng)};
    vertices.reserve(num_of_vertices);

    std::uniform_real_distribution<float> uniform{};

    for (auto i{0U}; i < num_of_vertices; ++i)
    {
        // transform Z to be the position later
        const auto vertex_position{uniform(rng) * (right() - left()) + left()};
        const auto mark{uniform(rng)};
        vertices.emplace_back(Point{mark, vertex_position});
    }
    return vertices;
}
