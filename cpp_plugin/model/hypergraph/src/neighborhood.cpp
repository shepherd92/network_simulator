#include <cassert>
#include <cmath>

#include "hypergraph_model.hpp"
#include "neighborhood.hpp"
#include "tools.hpp"

NeighborhoodPart::NeighborhoodPart(const float left, const float right)
    : left_{left}, right_{right}
{
}

Position NeighborhoodPart::left() const
{
    return left_;
}

Position NeighborhoodPart::right() const
{
    return right_;
}

void NeighborhoodPart::set_left(const Position value)
{
    left_ = value;
}

void NeighborhoodPart::set_right(const Position value)
{
    right_ = value;
}

Hyperbola::Hyperbola(
    const Position left,
    const Position right,
    const Position position,
    const Mark transformed_mark)
    : NeighborhoodPart{left, right}, position_{position}, transformed_mark_{transformed_mark}
{
}

std::optional<Hyperbola> Hyperbola::intersect_domain(const Position min_, const Position max_) const
{
    if (left() > max_ || right() < min_)
    {
        // no intersection
        return std::nullopt;
    }

    return Hyperbola(std::max(min_, left()), std::min(max_, right()), position(), transformed_mark());
}

std::vector<Hyperbola> Hyperbola::get_dominating_hyperbola_parts(const Hyperbola &other) const
{
    // assumption: domain of other contains the domain of this hyperbola
    std::vector<Hyperbola> dominating_hyperbolas{};

    if (less_than(other, left()) && less_than(other, right()))
    {
        // this hyperbola is less than the other in the entire domain
        dominating_hyperbolas.emplace_back(Hyperbola(left(), right(), other.position(), other.transformed_mark()));
        return dominating_hyperbolas;
    }
    if (!less_than(other, left()) && !less_than(other, right()))
    {
        // this hyperbola is greater than the other in the entire domain
        dominating_hyperbolas.emplace_back(*this);
        return dominating_hyperbolas;
    }

    // at this point we know that the hyperbolas intersect in the domain
    const auto V1{transformed_mark()};
    const auto V2{other.transformed_mark()};
    const auto y1{position()};
    const auto y2{other.position()};
    const bool hyperbolas_are_same_sided{(is_left_tail() && other.is_left_tail()) || (!is_left_tail() && !other.is_left_tail())};

    if (hyperbolas_are_same_sided)
    {
        const Position intersection{(V2 * y1 - V1 * y2) / (V2 - V1)};
        // larger transformed mark hyperbola dominates further away from the intersection
        if ((is_left_tail() && V1 > V2) || (!is_left_tail() && V1 < V2))
        {
            dominating_hyperbolas.emplace_back(Hyperbola(left(), intersection, y1, V1));
            dominating_hyperbolas.emplace_back(Hyperbola(intersection, right(), y2, V2));
        }
        else
        {
            dominating_hyperbolas.emplace_back(Hyperbola(left(), intersection, y2, V2));
            dominating_hyperbolas.emplace_back(Hyperbola(intersection, right(), y1, V1));
        }
    }
    else
    {
        const Position intersection{(V2 * y1 + V1 * y2) / (V2 + V1)};
        // right tail hyperbola dominates left from the intersection
        if (is_left_tail())
        {
            dominating_hyperbolas.emplace_back(Hyperbola(left(), intersection, y2, V2));
            dominating_hyperbolas.emplace_back(Hyperbola(intersection, right(), y1, V1));
        }
        else
        {
            dominating_hyperbolas.emplace_back(Hyperbola(left(), intersection, y1, V1));
            dominating_hyperbolas.emplace_back(Hyperbola(intersection, right(), y2, V2));
        }
    }

    return dominating_hyperbolas;
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
    return std::pow(std::abs(position() - x) / transformed_mark(), -1.F / gamma);
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

    const auto uniform_lower_limit{std::pow(std::abs(right() - position()), 1.F - 1.F / g)};
    const auto uniform_upper_limit{std::pow(std::abs(left() - position()), 1.F - 1.F / g)};
    std::uniform_real_distribution<float> uniform{uniform_lower_limit, uniform_upper_limit};

    const auto transformation_coefficient{left() < position() ? -1.F : 1.F};
    for (auto i{0U}; i < num_of_vertices; ++i)
    {
        // transform Z to be the position later
        const auto Z{uniform(rng)};
        const auto vertex_position{position() + transformation_coefficient * std::pow(Z, uniform_exponent)};
        // avoid marks larger than 1 due to floating point errors
        const auto max_mark{std::min((*this)(vertex_position, g), 1.F)};
        const auto mark{max_mark * uniform_01(rng)};
        vertices.emplace_back(Point{mark, vertex_position});
    }
    return vertices;
}

float Hyperbola::integral(const float gamma) const
{
    return gamma / (1. - gamma) *
           std::pow(transformed_mark(), 1. / gamma) *
           std::abs((std::pow(std::abs(position() - right()), -(1. - gamma) / gamma) -
                     std::pow(std::abs(position() - left()), -(1. - gamma) / gamma)));
}

bool Hyperbola::less_than(const Hyperbola &other, const Position x) const
{
    return std::abs(x - position()) / transformed_mark() > std::abs(x - other.position()) / other.transformed_mark();
}

bool Hyperbola::is_left_tail() const
{
    return right() < position();
}

Center::Center(const Position left, const Position right)
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

Position Hyperbola::position() const
{
    return position_;
}

Mark Hyperbola::transformed_mark() const
{
    return transformed_mark_;
}
