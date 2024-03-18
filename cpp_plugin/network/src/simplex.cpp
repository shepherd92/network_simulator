#include <cassert>
#include <iostream>
#include <map>
#include <mutex>

#include "simplex.h"
#include "tools.h"

std::size_t SimplexHash::operator()(const Simplex &simplex) const
{
    const auto &vertices{simplex.vertices()};
    std::size_t seed{vertices.size()};
    std::for_each(
        std::execution::seq,
        vertices.begin(), vertices.end(),
        [&seed](const auto &vertex)
        {
            seed ^= vertex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        });

    return seed;
}

Simplex::Simplex(const PointIdList &vertices)
    : vertices_{vertices}
{
    assert(!vertices_.empty());
    std::sort(vertices_.begin(), vertices_.end());
    for (const auto &vertex_id : vertices_)
    {
        bloom_filter_.set(vertex_id & 63);
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

SimplexList Simplex::skeleton(const Dimension max_dimension) const
{
    SimplexList result{};
    for (Dimension dimension{0}; dimension <= max_dimension; ++dimension)
    {
        const auto faces_{faces(dimension)};
        std::copy(faces_.begin(), faces_.end(), std::back_inserter(result));
    }
    return result;
}

SimplexList Simplex::faces(const Dimension max_dimension) const
{
    if (dimension() <= max_dimension)
    {
        return SimplexList{};
    }

    PointIdList current_combination(max_dimension + 1U);
    SimplexList result{};
    combination_util(max_dimension, 0U, result, current_combination, 0U);
    return result;
}

void Simplex::combination_util(
    const Dimension max_dimension,
    const uint32_t combination_index,
    SimplexList &result,
    PointIdList &current_combination,
    const uint32_t array_index) const
{
    if (combination_index == max_dimension + 1U)
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
    combination_util(max_dimension, combination_index + 1U, result, current_combination, array_index + 1U);
    combination_util(max_dimension, combination_index, result, current_combination, array_index + 1U);
}

SimplexList create_simplices(const ISimplexList &simplices_in)
{
    SimplexList simplices;
    std::mutex mutex{};
    const auto total{simplices_in.size()};
    std::atomic<uint32_t> counter{0U};
    simplices.reserve(total);

    std::for_each(
        execution_policy,
        simplices_in.begin(),
        simplices_in.end(),
        [&](const auto &simplex_in)
        {
            std::lock_guard<std::mutex> lock(mutex);
            simplices.emplace_back(Simplex(simplex_in));
        });

    return simplices;
}

ISimplexList create_raw_simplices(const SimplexList &simplices_in)
{
    ISimplexList simplices;
    for (const auto &simplex : simplices_in)
    {
        simplices.push_back(simplex.vertices());
    }
    return simplices;
}

SimplexList get_skeleton_simplices(const SimplexList &simplices, const Dimension max_dimension)
{
    SimplexList skeleton_simplices{};
    std::mutex mutex{};

    for (Dimension dimension{0}; dimension <= max_dimension; ++dimension)
    {
        const auto faces{get_faces_simplices(simplices, dimension)};
        std::copy(faces.begin(), faces.end(), std::back_inserter(skeleton_simplices));
    }
    return skeleton_simplices;
}

SimplexList get_faces_simplices(const SimplexList &simplices, const Dimension dimension)
{
    SimplexSet faces_set{};
    std::mutex mutex{};

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            const auto faces{simplex.faces(dimension)};
            std::lock_guard<std::mutex> lock(mutex);
            std::copy(faces.begin(), faces.end(), std::inserter(faces_set, faces_set.end()));
        });

    SimplexList faces{};
    faces.reserve(faces_set.size());
    for (auto it = faces_set.begin(); it != faces_set.end();)
    {
        faces.push_back(std::move(faces_set.extract(it++).value()));
    }

    sort_simplices(faces, true);
    return faces;
}

SimplexList get_cofaces_simplices(const SimplexList &simplices_in, const Simplex &simplex)
{
    SimplexList cofaces{};

    std::copy_if(
        execution_policy,
        simplices_in.begin(),
        simplices_in.end(),
        std::back_inserter(cofaces),
        [&](const auto &other_simplex)
        { return simplex.is_face(other_simplex); });

    return cofaces;
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

SimplexList select_higher_dimensional_simplices(const SimplexList &simplices, const Dimension min_dimension)
{
    SimplexList selected_simplices;

    std::copy_if(simplices.begin(), simplices.end(),
                 std::back_inserter(selected_simplices),
                 [min_dimension](const Simplex &simplex)
                 { return simplex.dimension() >= min_dimension; });

    sort_simplices(selected_simplices, true);
    return selected_simplices;
}

void sort_simplices(SimplexList &simplices, const bool ascending)
{
    if (ascending)
    {
        std::sort(simplices.begin(), simplices.end(),
                  [](const Simplex &lhs, const Simplex &rhs)
                  { return lhs.vertices().size() < rhs.vertices().size(); });
    }
    else
    {
        std::sort(simplices.begin(), simplices.end(),
                  [](const Simplex &lhs, const Simplex &rhs)
                  { return lhs.vertices().size() > rhs.vertices().size(); });
    }
}

std::vector<Dimension> calc_dimension_distribution(const ISimplexList &simplices_in)
{
    const auto simplices{create_simplices(simplices_in)};
    return calc_dimension_distribution(simplices);
}

SimplexList filter_simplices(const SimplexList &simplices, const PointIdList &vertices_to_keep)
{
    SimplexList filtered_simplices{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            auto keep_simplex{true};
            for (auto vertex : simplex.vertices())
            {
                const auto vertex_is_kept{std::find(vertices_to_keep.begin(), vertices_to_keep.end(), vertex) != vertices_to_keep.end()};
                if (!vertex_is_kept)
                {
                    keep_simplex = false;
                    break;
                }
            }
            if (keep_simplex)
            {
                std::lock_guard<std::mutex> lock{mutex};
                filtered_simplices.push_back(simplex);
            }
        });
    return filtered_simplices;
}

std::vector<Dimension> calc_dimension_distribution(const SimplexList &simplices)
{
    std::mutex mutex{};
    const auto total{simplices.size()};
    std::atomic<uint32_t> counter{0U};
    std::vector<Dimension> result{};
    result.reserve(total);

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            const auto dimension{simplex.dimension()};
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back(dimension);
            log_progress(++counter, total, 1000U, "Calc simplex dimension distribution");
        });
    log_progress(counter, total, 1U, "Calc simplex dimension distribution");

    std::sort(result.begin(), result.end());
    return result;
}