#include <cassert>
#include <map>
#include <vector>

#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

SimplexList::SimplexList(const ISimplexList &input)
    : simplices_{}
{
    simplices_.reserve(input.size());
    std::mutex mutex{};
    const auto total{input.size()};
    std::atomic<uint32_t> counter{0U};
    simplices_.reserve(total);

    std::for_each(
        execution_policy,
        input.begin(),
        input.end(),
        [&](const auto &simplex_in)
        {
            std::lock_guard<std::mutex> lock(mutex);
            simplices_.emplace_back(Simplex(simplex_in));
        });
    sort(true);
}

SimplexList::SimplexList(const std::vector<Simplex> &input)
    : simplices_{input}
{
    sort(true);
}

SimplexList::SimplexList(SimplexSet &input)
    : simplices_{}
{
    simplices_.reserve(input.size());
    for (auto it = input.begin(); it != input.end();)
    {
        simplices_.push_back(std::move(input.extract(it++).value()));
    }
    sort(true);
}

ISimplexList SimplexList::raw() const
{
    ISimplexList result;
    for (const auto &simplex : simplices_)
    {
        result.push_back(simplex.vertices());
    }
    return result;
}

std::vector<Dimension> SimplexList::calc_dimension_distribution() const
{
    std::mutex mutex{};
    const auto total{simplices_.size()};
    std::atomic<uint32_t> counter{0U};
    std::vector<Dimension> result{};
    result.reserve(total);

    std::for_each(
        execution_policy,
        simplices_.begin(),
        simplices_.end(),
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

SimplexList SimplexList::filter(const PointIdList &vertices_to_keep) const
{
    std::vector<Simplex> filtered_simplices{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices_.begin(),
        simplices_.end(),
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
    return SimplexList{filtered_simplices};
}

void SimplexList::sort(const bool ascending)
{
    if (ascending)
    {
        std::sort(simplices_.begin(), simplices_.end(),
                  [](const Simplex &lhs, const Simplex &rhs)
                  { return lhs.vertices().size() < rhs.vertices().size(); });
    }
    else
    {
        std::sort(simplices_.begin(), simplices_.end(),
                  [](const Simplex &lhs, const Simplex &rhs)
                  { return lhs.vertices().size() > rhs.vertices().size(); });
    }
}

SimplexList SimplexList::select_higher_dimensional_simplices(const Dimension min_dimension) const
{
    std::vector<Simplex> selected_simplices{};
    std::copy_if(simplices_.begin(), simplices_.end(),
                 std::back_inserter(selected_simplices),
                 [min_dimension](const Simplex &simplex)
                 { return simplex.dimension() >= min_dimension; });

    return SimplexList{selected_simplices};
}

SimplexList SimplexList::faces(const Dimension dimension) const
{
    SimplexSet faces_set{};
    std::mutex mutex{};

    std::for_each(
        execution_policy,
        simplices_.begin(),
        simplices_.end(),
        [&](const auto &simplex)
        {
            const auto faces{simplex.faces(dimension)};
            std::lock_guard<std::mutex> lock(mutex);
            std::copy(faces.begin(), faces.end(), std::inserter(faces_set, faces_set.end()));
        });

    return SimplexList{faces_set};
}

SimplexList SimplexList::cofaces(const Simplex &simplex) const
{
    std::vector<Simplex> cofaces{};

    std::copy_if(
        execution_policy,
        simplices_.begin(),
        simplices_.end(),
        std::back_inserter(cofaces),
        [&](const auto &other_simplex)
        { return simplex.is_face(other_simplex); });

    return SimplexList{cofaces};
}

SimplexList SimplexList::facets() const
{
    std::mutex mutex{};
    std::vector<Simplex> facets{};
    const auto total{simplices_.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        simplices_.begin(),
        simplices_.end(),
        [&](auto &&simplex)
        {
            const auto first_it{simplices_.begin() + (&simplex - &(simplices_[0]))};
            auto first_is_face{false};
            for (auto second_it{first_it + 1}; second_it < simplices_.end(); ++second_it)
            {
                first_is_face = false;
                if (simplex.is_face(*second_it))
                {
                    first_is_face = true;
                    break;
                }
            }

            if (!first_is_face)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                facets.push_back(simplex);
            }
            log_progress(++counter, total, 1000U, "Calc facets");
        });

    log_progress(counter, total, 1U, "Calc facets");
    return SimplexList{facets};
}

SimplexList SimplexList::simplices_by_dimension(const Dimension dimension) const
{
    std::vector<Simplex> selected_simplices;

    std::copy_if(simplices_.begin(), simplices_.end(),
                 std::back_inserter(selected_simplices),
                 [dimension](const Simplex &simplex)
                 { return simplex.dimension() == dimension; });

    return SimplexList{selected_simplices};
}

SimplexList SimplexList::skeleton(const Dimension max_dimension) const
{
    SimplexList skeleton_simplices{};
    std::mutex mutex{};

    for (Dimension dimension{0}; dimension <= max_dimension; ++dimension)
    {
        const auto faces_dim{faces(dimension)};
        skeleton_simplices += faces_dim;
    }
    return skeleton_simplices;
}

uint32_t SimplexList::size() const
{
    return simplices_.size();
}

SimplexList SimplexList::operator+(const SimplexList &other) const
{
    SimplexList result{*this};
    std::copy(other.simplices_.begin(), other.simplices_.end(), std::back_inserter(result.simplices_));
    return result;
}

void SimplexList::operator+=(const SimplexList &other)
{
    std::copy(other.simplices_.begin(), other.simplices_.end(), std::back_inserter(simplices_));
}

std::vector<PointId> SimplexList::vertices() const
{
    std::unordered_set<PointId> unique_vertices{};
    std::for_each(
        std::execution::seq,
        simplices_.begin(),
        simplices_.end(),
        [&](const auto &interaction)
        {
            for (const auto &vertex : interaction.vertices())
            {
                unique_vertices.insert(vertex);
            }
        });
    PointIdList vertices{unique_vertices.begin(), unique_vertices.end()};
    std::sort(vertices.begin(), vertices.end());
    return vertices;
}

const std::vector<Simplex> &SimplexList::simplices() const
{
    return simplices_;
}

std::map<PointId, std::vector<int32_t>> SimplexList::vertex_simplex_map() const
{
    std::map<PointId, std::vector<int32_t>> vertex_simplex_map{};
    for (const auto &vertex : vertices())
    {
        vertex_simplex_map[vertex] = std::vector<int32_t>{};
    }

    for (auto index{0U}; index < simplices_.size(); ++index)
    {
        for (const auto &vertex : simplices_[index].vertices())
        {
            vertex_simplex_map[vertex].push_back(index);
        }
    }
    // note: assume that the vertex_simplex_map is sorted
    return vertex_simplex_map;
}

std::vector<uint32_t> SimplexList::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension) const
{
    assert(neighbor_dimension > simplex_dimension);
    const auto &simplices{simplices_by_dimension(simplex_dimension).simplices()};
    const auto &possible_cofaces{simplices_by_dimension(neighbor_dimension).simplices()};
    std::vector<uint32_t> degree_sequence{};
    degree_sequence.reserve(simplices.size());

    std::mutex mutex{};
    const auto total{simplices.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        std::execution::seq,
        simplices.begin(),
        simplices.end(),
        [&](auto &&simplex)
        {
            std::atomic<uint32_t> degree{0U};

            std::for_each(
                execution_policy,
                possible_cofaces.begin(),
                possible_cofaces.end(),
                [&](const auto &neighbor)
                {
                    if (simplex.is_face(neighbor))
                    {
                        ++degree;
                    }
                });

            std::lock_guard<std::mutex> lock_guard(mutex);
            degree_sequence.push_back(degree);
            log_progress(++counter, total, 1000U, "Calc degree sequence");
        });
    log_progress(counter, total, 1U, "Calc degree sequence");

    return degree_sequence;
}

SimplexList::iterator SimplexList::begin()
{
    return simplices_.begin();
}

SimplexList::const_iterator SimplexList::begin() const
{
    return simplices_.begin();
}

SimplexList::const_iterator SimplexList::cbegin() const
{
    return simplices_.cbegin();
}

SimplexList::iterator SimplexList::end()
{
    return simplices_.end();
}

SimplexList::const_iterator SimplexList::end() const
{
    return simplices_.end();
}

SimplexList::const_iterator SimplexList::cend() const
{
    return simplices_.cend();
}
