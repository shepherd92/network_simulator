#include <atomic>
#include <execution>
#include <mutex>

#include "network.h"
#include "numpy_cpp_conversion.h"

Network::Network(const Dimension max_dimension)
    : max_dimension_{max_dimension},
      simplex_tree_{},
      interactions_{std::nullopt},
      facets_{std::nullopt},
      persistent_cohomology_{simplex_tree_},
      is_persistence_calculated_{false}
{
}

void Network::add_simplices(const ISimplexList &simplices)
{
    for (const auto &simplex : simplices)
    {
        add_simplex(simplex);
    }
    is_persistence_calculated_ = false;
}

void Network::reset()
{
    simplex_tree_ = SimplexTree{};
    is_persistence_calculated_ = false;
    interactions_ = std::nullopt;
    facets_ = std::nullopt;
}

void Network::add_simplex(const ISimplex &simplex)
{
    if (simplex.size() - 1 <= max_dimension_)
    {
        simplex_tree_.insert_simplex_and_subfaces(simplex);
    }
    else
    {
        ISimplex current_combination(max_dimension_ + 1U);
        ISimplexList result{};
        combination_util(simplex, 0U, result, current_combination, 0U);
        for (const auto &skeleton_simplex : result)
        {
            simplex_tree_.insert_simplex_and_subfaces(skeleton_simplex);
        }
    }
}

void Network::combination_util(
    const ISimplex &simplex,
    const uint32_t combination_index,
    ISimplexList &result,
    ISimplex &current_combination,
    const uint32_t array_index) const
{
    if (combination_index == max_dimension_ + 1U)
    {
        // combination ready
        result.push_back(current_combination);
        return;
    }

    if (array_index > simplex.size() - 1U)
        return;

    current_combination[combination_index] = simplex[array_index];
    combination_util(simplex, combination_index + 1U, result, current_combination, array_index + 1U);
    combination_util(simplex, combination_index, result, current_combination, array_index + 1U);
}

PersistentCohomology &Network::get_persistent_cohomology()
{
    if (!is_persistence_calculated_)
    {
        calc_persistent_cohomology();
        is_persistence_calculated_ = true;
    }
    return persistent_cohomology_;
}

void Network::calc_persistent_cohomology()
{
    persistent_cohomology_.init_coefficients(2);
    std::cout << __LINE__ << std::endl;
    persistent_cohomology_.compute_persistent_cohomology();
    std::cout << __LINE__ << std::endl;
}

std::vector<ISimplexList> Network::calc_persistence_pairs()
{
    const auto &persistent_cohomology{get_persistent_cohomology()};
    const auto &persistent_pairs{persistent_cohomology.get_persistent_pairs()};
    std::vector<ISimplexList> result{};
    result.reserve(persistent_pairs.size());
    for (const auto &persistent_interval : persistent_pairs)
    {
        ISimplex birth_simplex{};
        for (auto vertex : simplex_tree_.simplex_vertex_range(std::get<0>(persistent_interval)))
        {
            birth_simplex.push_back(vertex);
        }
        ISimplex death_simplex{};
        for (auto vertex : simplex_tree_.simplex_vertex_range(std::get<1>(persistent_interval)))
        {
            death_simplex.push_back(vertex);
        }
        result.push_back({birth_simplex, death_simplex});
    }

    return result;
}

std::vector<int32_t> Network::calc_betti_numbers()
{
    const auto &persistent_cohomology{get_persistent_cohomology()};
    const auto result{persistent_cohomology.betti_numbers()};
    return result;
}

void Network::expand()
{
    simplex_tree_.expansion(max_dimension_ + 1U);
    is_persistence_calculated_ = false;
}

uint32_t Network::num_simplices()
{
    return simplex_tree_.num_simplices();
}

uint32_t Network::num_vertices() const
{
    return simplex_tree_.num_vertices();
}

ISimplexList Network::get_facets()
{
    ISimplexList result{};
    if (!facets_)
    {
        facets_ = calc_facets();
    }
    return *facets_;
}

ISimplexList Network::calc_facets()
{
    auto simplices{get_simplices()};
    std::mutex mutex;
    ISimplexList facets{};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        std::execution::seq,
        simplices.begin(),
        simplices.end(),
        [&](auto &&simplex)
        {
            const auto first_it{simplices.begin() + (&simplex - &simplices[0])};
            if (++counter % 10000 == 0 && simplices.size() > 100000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rC++: Extracting facets ... " << counter << " / " << simplices.size();
            }
            auto first_is_face{false};
            for (auto second_it{first_it + 1}; second_it < simplices.end(); ++second_it)
            {
                first_is_face = false;
                if (std::includes(
                        second_it->begin(), second_it->end(),
                        first_it->begin(), first_it->end()))
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
        });

    return facets;
}

void Network::set_facets(const ISimplexList &facets)
{
    facets_ = facets;
}

ISimplexList Network::get_simplices()
{
    return convert_to_raw_simplices<>(simplex_tree_.filtration_simplex_range());
}

ISimplexList Network::get_skeleton(const Dimension max_dimension)
{
    return convert_to_raw_simplices<>(simplex_tree_.skeleton_simplex_range(max_dimension));
}

template <typename Iterator>
ISimplexList Network::convert_to_raw_simplices(const Iterator &simplex_range)
{
    ISimplexList result{};
    for (const auto &simplex : simplex_range)
    {
        ISimplex vertices{};
        vertices.reserve(simplex_tree_.dimension(simplex) + 1U);
        for (auto vertex : simplex_tree_.simplex_vertex_range(simplex))
        {
            vertices.push_back(vertex);
        }
        std::sort(vertices.begin(), vertices.end());
        result.push_back(vertices);
    }
    return result;
}

void Network::set_simplices(const ISimplexList &simplices)
{
    simplex_tree_ = SimplexTree();
    add_simplices(simplices);
    is_persistence_calculated_ = false;
}

void Network::keep_only_vertices(const VertexList &vertices)
{
    filter_simplex_tree(vertices);
    filter_interactions(vertices);
    facets_ = std::nullopt;
    is_persistence_calculated_ = false;
}

std::vector<uint32_t> Network::calc_simplex_dimension_distribution()
{
    std::vector<uint32_t> result(max_dimension_ + 1U, 0U);

    for (const auto &simplex : simplex_tree_.filtration_simplex_range())
    {
        ++result[simplex_tree_.dimension(simplex)];
    }
    return result;
}

std::vector<uint32_t> Network::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    std::vector<uint32_t> result{};
    for (const auto &simplex : simplex_tree_.filtration_simplex_range())
    {
        if (simplex_tree_.dimension(simplex) == static_cast<int32_t>(simplex_dimension))
            result.push_back(
                simplex_tree_.cofaces_simplex_range(simplex, neighbor_dimension - simplex_dimension).size());
    }
    std::sort(result.begin(), result.end());
    return result;
}

void Network::filter_simplex_tree(const VertexList &vertices)
{
    SimplexTree filtered_simplex_tree{};
    for (const auto &simplex : simplex_tree_.filtration_simplex_range())
    {
        auto keep_simplex{true};
        for (auto vertex : simplex_tree_.simplex_vertex_range(simplex))
        {
            const auto vertex_in_simplex{std::find(vertices.begin(), vertices.end(), vertex) != vertices.end()};
            if (!vertex_in_simplex)
            {
                keep_simplex = false;
                break;
            }
        }
        if (keep_simplex)
        {

            filtered_simplex_tree.insert_simplex_and_subfaces(simplex_tree_.simplex_vertex_range(simplex));
        }
    }
    simplex_tree_ = filtered_simplex_tree;
}

void Network::filter_interactions(const ISimplex &vertices)
{
    (void)vertices;
}

Dimension Network::get_max_dimension() const
{
    return max_dimension_;
}

void Network::set_max_dimension(const Dimension max_dimension)
{
    max_dimension_ = max_dimension;
}
