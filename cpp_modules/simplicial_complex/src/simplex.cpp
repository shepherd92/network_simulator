#include "simplex.h"

Simplex::Simplex(const std::vector<int32_t> &vertices)
    : vertices_(std::set<int32_t>(vertices.begin(), vertices.end()))
{
}

const std::set<int32_t> &Simplex::vertices() const
{
    return vertices_;
}

std::vector<int32_t> Simplex::operator-(const Simplex &other) const
{
    std::vector<int32_t> remaining_vertices;
    std::set_difference(
        vertices_.begin(), vertices_.end(),
        other.vertices().begin(), other.vertices().end(),
        std::inserter(remaining_vertices, remaining_vertices.begin()));
    return remaining_vertices;
}

uint32_t Simplex::dimension() const
{
    return vertices_.size() - 1;
}

bool Simplex::is_face(const Simplex &other) const
{
    return std::includes(other.vertices().begin(), other.vertices().end(), vertices_.begin(), vertices_.end());
}

std::vector<Simplex> create_simplices(const std::vector<std::vector<int32_t>> &simplices_in)
{
    std::vector<Simplex> simplices;
    for (const auto &simplex : simplices_in)
    {
        simplices.push_back(Simplex(simplex));
    }
    return simplices;
}

std::vector<Simplex> select_simplices_by_dimension(
    const std::vector<Simplex> &simplices,
    const uint32_t dimension)
{
    std::vector<Simplex> selected_simplices;

    std::copy_if(simplices.begin(), simplices.end(),
                 std::back_inserter(selected_simplices),
                 [dimension](const Simplex &simplex)
                 { return simplex.dimension() == dimension; });

    return selected_simplices;
}

std::vector<Simplex> select_higher_dimensional_simplices(
    const std::vector<Simplex> &simplices,
    const uint32_t dimension)
{
    std::vector<Simplex> selected_simplices;

    std::copy_if(simplices.begin(), simplices.end(),
                 std::back_inserter(selected_simplices),
                 [dimension](const Simplex &simplex)
                 { return simplex.dimension() >= dimension; });

    return selected_simplices;
}

std::vector<Simplex> sort_simplices(
    const std::vector<Simplex> &simplices,
    const bool ascending)
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
