#include <iostream>

#include "simplex.h"

std::size_t Hash::operator()(const Simplex &simplex) const
{
    const auto &vertices = simplex.vertices();
    std::size_t seed{vertices.size()};
    for (auto &vertex : vertices)
    {
        seed ^= vertex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

Simplex::Simplex(const VertexList &vertices)
    : vertices_(vertices)
{
}

const VertexList &Simplex::vertices() const
{
    return vertices_;
}

Simplex Simplex::operator-(const Simplex &other) const
{
    VertexList remaining_vertices;
    std::set_difference(
        vertices_.begin(), vertices_.end(),
        other.vertices().begin(), other.vertices().end(),
        std::inserter(remaining_vertices, remaining_vertices.begin()));
    return Simplex(remaining_vertices);
}

bool Simplex::operator==(const Simplex &other) const
{
    return std::is_permutation(vertices_.begin(), vertices_.end(), other.vertices().begin(), other.vertices().end());
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

bool Simplex::is_face(const Simplex &other) const
{
    return std::includes(other.vertices().begin(), other.vertices().end(), vertices_.begin(), vertices_.end());
}

SimplexList Simplex::get_skeleton(const Dimension max_dimension) const
{
    if (dimension() <= max_dimension)
    {
        return SimplexList{*this};
    }
    VertexList current_combination(max_dimension + 1U);
    SimplexList result{};

    combination_util(max_dimension, 0U, result, current_combination, 0U);
    return result;
}

void Simplex::combination_util(
    const Dimension max_dimension,
    const uint32_t combination_index,
    SimplexList &result,
    VertexList &current_combination,
    const uint32_t array_index) const
{
    if (combination_index == max_dimension + 1U)
    {
        // combination ready
        result.push_back(Simplex{current_combination});
        return;
    }

    if (array_index > dimension())
        return;

    current_combination[combination_index] = vertices()[array_index];
    combination_util(max_dimension, combination_index + 1U, result, current_combination, array_index + 1U);
    combination_util(max_dimension, combination_index, result, current_combination, array_index + 1U);
}

SimplexList create_simplices(const std::vector<VertexList> &simplices_in)
{
    SimplexList simplices;
    for (const auto &simplex : simplices_in)
    {
        simplices.push_back(Simplex(simplex));
    }
    return simplices;
}

SimplexList select_simplices_by_dimension(const SimplexList &simplices, const Dimension dimension)
{
    SimplexList selected_simplices;

    std::copy_if(simplices.begin(), simplices.end(),
                 std::back_inserter(selected_simplices),
                 [dimension](const Simplex &simplex)
                 { return simplex.dimension() == dimension; });

    return selected_simplices;
}

SimplexList select_higher_dimensional_simplices(const SimplexList &simplices, const Dimension dimension)
{
    SimplexList selected_simplices;

    std::copy_if(simplices.begin(), simplices.end(),
                 std::back_inserter(selected_simplices),
                 [dimension](const Simplex &simplex)
                 { return simplex.dimension() >= dimension; });

    return selected_simplices;
}

SimplexList sort_simplices(const SimplexList &simplices, const bool ascending)
{
    auto sorted_simplices{simplices};
    if (ascending)
    {
        std::sort(sorted_simplices.begin(), sorted_simplices.end(),
                  [](const Simplex &lhs, const Simplex &rhs)
                  { return lhs.vertices().size() < rhs.vertices().size(); });
    }
    else
    {
        std::sort(sorted_simplices.begin(), sorted_simplices.end(),
                  [](const Simplex &lhs, const Simplex &rhs)
                  { return lhs.vertices().size() > rhs.vertices().size(); });
    }
    return sorted_simplices;
}