#include <iostream>
#include <map>
#include <mutex>

#include "combinations.h"
#include "simplex.h"
#include "tools.h"

void Simplex::combination_util(
    const Dimension dimension_,
    const uint32_t combination_index,
    std::vector<Simplex> &result,
    PointIdList &current_combination,
    const uint32_t array_index) const
{
    if (combination_index == dimension_ + 1U)
    {
        // combination ready
        result.emplace_back(Simplex{current_combination});
        return;
    }

    if (static_cast<int32_t>(array_index) > dimension())
    {
        return;
    }

    current_combination[combination_index] = vertices()[array_index];
    combination_util(dimension_, combination_index + 1U, result, current_combination, array_index + 1U);
    combination_util(dimension_, combination_index, result, current_combination, array_index + 1U);
}

uint64_t SimplexHash::operator()(const Simplex &simplex) const
{
    uint64_t seed{simplex.vertices().size()};
    for (auto vertex : simplex.vertices())
    {
        vertex = ((vertex >> 16) ^ vertex) * 0x45d9f3b;
        vertex = ((vertex >> 16) ^ vertex) * 0x45d9f3b;
        vertex = (vertex >> 16) ^ vertex;
        seed ^= vertex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

Simplex::Simplex(const PointIdList &vertices)
    : vertices_{vertices}, bloom_filter_{}
{
    vertices_.shrink_to_fit();
    std::sort(vertices_.begin(), vertices_.end());

    for (const auto vertex : vertices_)
    {
        bloom_filter_.set(vertex % BLOOM_FILTER_SIZE);
    }
}

const PointIdList &Simplex::vertices() const
{
    return vertices_;
}

Simplex Simplex::operator-(const Simplex &other) const
{
    PointIdList remaining_vertices;
    std::set_difference(
        vertices_.begin(), vertices_.end(),
        other.vertices().begin(), other.vertices().end(),
        std::inserter(remaining_vertices, remaining_vertices.begin()));
    return Simplex{remaining_vertices};
}

bool Simplex::operator==(const Simplex &other) const
{
    return vertices_ == other.vertices();
}

bool Simplex::operator<(const Simplex &other) const
{
    return vertices_ < other.vertices();
}

std::ostream &operator<<(std::ostream &os, const Simplex &simplex)
{
    for (const auto &vertex : simplex.vertices())
    {
        os << vertex << " ";
    }
    return os;
}

Dimension Simplex::dimension() const
{
    return vertices_.size() - 1;
}

std::vector<Simplex> Simplex::skeleton(const Dimension max_dimension) const
{
    std::vector<Simplex> result{};
    for (Dimension dimension{0}; dimension <= max_dimension; ++dimension)
    {
        const auto faces_{faces(dimension)};
        std::copy(faces_.begin(), faces_.end(), std::back_inserter(result));
    }
    return result;
}

std::vector<Simplex> Simplex::faces(const Dimension dimension_param) const
{
    std::vector<Simplex> result{};
    if (dimension_param > dimension())
    {
        return result;
    }

    result.reserve(binomial_coefficient(dimension() + 1, dimension_param + 1));
    auto vertex_ids{vertices()};

    for_each_combination(
        vertex_ids.begin(),
        vertex_ids.begin() + dimension_param + 1,
        vertex_ids.end(),
        [&result](PointIdList::iterator first, PointIdList::iterator last)
        {
            const PointIdList face_vertices{first, last};
            result.emplace_back(Simplex{face_vertices});
            return false;
        });

    return result;
}
