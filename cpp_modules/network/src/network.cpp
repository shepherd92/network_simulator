#include <atomic>
#include <execution>
#include <mutex>

#include "network.h"
#include "numpy_cpp_conversion.h"
#include "simplex.h"

Network::Network(const Dimension max_dimension)
    : max_dimension_{max_dimension},
      vertices_{},
      simplex_tree_{std::nullopt},
      interactions_{std::nullopt},
      facets_{std::nullopt},
      persistent_cohomology_{nullptr}
{
}

Network::~Network()
{
    reset_persistence();
}

void Network::reset()
{
    reset_simplex_tree();
    interactions_ = std::nullopt;
    facets_ = std::nullopt;
}

void Network::reset_simplex_tree()
{
    simplex_tree_ = std::nullopt;
    reset_persistence();
}

void Network::reset_persistence()
{
    delete persistent_cohomology_;
    persistent_cohomology_ = nullptr;
}

void Network::add_simplices_interface(const ISimplexList &simplices)
{
    std::mutex mutex;
    ISimplexList skeleton_simplices{get_skeleton(simplices)};

    if (!simplex_tree_.has_value())
    {
        simplex_tree_ = SimplexTree{};
    }

    std::for_each(
        std::execution::seq,
        skeleton_simplices.begin(),
        skeleton_simplices.end(),
        [this](const auto &simplex)
        {
            simplex_tree_->insert_simplex_and_subfaces(simplex);
        });
    reset_persistence();
}

void Network::add_simplices(const SimplexList &simplices)
{
    SimplexList skeleton_simplices{get_skeleton_simplices(simplices, max_dimension_)};

    if (!simplex_tree_.has_value())
    {
        simplex_tree_ = SimplexTree{};
    }

    std::for_each(
        std::execution::seq,
        skeleton_simplices.begin(),
        skeleton_simplices.end(),
        [this](const auto &simplex)
        {
            simplex_tree_->insert_simplex_and_subfaces(simplex.vertices());
        });
    reset_persistence();
}

void Network::add_vertices(const VertexList &vertices)
{
    if (!simplex_tree_.has_value())
    {
        simplex_tree_ = SimplexTree{};
    }
    std::for_each(
        std::execution::seq,
        vertices.begin(),
        vertices.end(),
        [this](const auto &vertex)
        {
            simplex_tree_->insert_simplex({vertex});
        });
    reset_persistence();
}

const PersistentCohomology &Network::get_persistence()
{
    if (!persistent_cohomology_)
    {
        calc_persistent_cohomology();
    }
    return *persistent_cohomology_;
}

void Network::calc_persistent_cohomology()
{
    auto &simplex_tree{get_simplex_tree()};
    reset_persistence();
    persistent_cohomology_ = new PersistentCohomology{simplex_tree};
    persistent_cohomology_->init_coefficients(2);
    persistent_cohomology_->compute_persistent_cohomology();
}

SimplexTree &Network::get_simplex_tree()
{
    assert(simplex_tree_.has_value() || interactions_.has_value());
    if (!simplex_tree_.has_value())
    {
        create_simplex_tree_from_interactions();
    }
    return *simplex_tree_;
}

const VertexList &Network::get_vertices_interface() const
{
    return vertices_;
}

void Network::set_vertices(const VertexList &vertices)
{
    vertices_ = vertices;
}

VertexList Network::get_vertices(const SimplexHandle &simplex_handle)
{
    assert(simplex_tree_.has_value());
    VertexList result{};
    if (simplex_handle != simplex_tree_->null_simplex())
    {
        for (const auto &vertex : simplex_tree_->simplex_vertex_range(simplex_handle))
        {
            result.push_back(vertex);
        }
    }
    return result;
}

std::vector<ISimplexList> Network::calc_persistence_pairs()
{
    const auto &persistent_cohomology{get_persistence()};
    const auto &persistent_pairs{persistent_cohomology.get_persistent_pairs()};
    std::vector<ISimplexList> result{};
    result.reserve(persistent_pairs.size());
    std::mutex mutex{};

    std::for_each(
        execution_policy,
        persistent_pairs.begin(),
        persistent_pairs.end(),
        [this, &result, &mutex](const auto &persistent_interval)
        {
            VertexList birth_simplex{get_vertices(std::get<0>(persistent_interval))};
            VertexList death_simplex{get_vertices(std::get<1>(persistent_interval))};
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back({birth_simplex, death_simplex});
        });

    return result;
}

std::vector<int32_t> Network::calc_betti_numbers()
{
    const auto &persistent_cohomology{get_persistence()};
    const auto result{persistent_cohomology.betti_numbers()};
    return result;
}

void Network::expand()
{
    assert(simplex_tree_.has_value());
    simplex_tree_->expansion(max_dimension_ + 1U);
    reset_persistence();
}

ISimplexList Network::get_simplices_interface()
{
    return convert_to_raw_simplices<>(get_simplices());
}

ISimplexList Network::get_interactions_interface() const
{
    return create_raw_simplices(get_interactions());
}

ISimplexList Network::get_facets_interface()
{
    return create_raw_simplices(get_facets());
}

uint32_t Network::num_vertices()
{
    return vertices_.size();
}

std::vector<Dimension> Network::calc_simplex_dimension_distribution()
{
    const auto &simplices{get_simplices()};
    const auto result{calc_dimension_distribution(simplices)};
    return result;
}

std::vector<Dimension> Network::calc_interaction_dimension_distribution()
{
    assert(interactions_ && "No interactions");
    const auto result{calc_dimension_distribution(*interactions_)};
    return result;
}

std::vector<Dimension> Network::calc_facet_dimension_distribution()
{
    const auto &facets{get_facets()};
    const auto result{calc_dimension_distribution(facets)};
    return result;
}

std::vector<Dimension> Network::calc_dimension_distribution(const SimplexList &simplices) const
{
    std::vector<Dimension> result{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back(simplex.dimension());
        });
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<Dimension> Network::calc_dimension_distribution(const ISimplexList &simplices) const
{
    std::vector<Dimension> result{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back(simplex.size() - 1U);
        });
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<Dimension> Network::calc_dimension_distribution(const SimplexHandleList &simplices)
{
    auto &simplex_tree(get_simplex_tree());
    std::vector<Dimension> result{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back(simplex_tree.dimension(simplex));
        });
    std::sort(result.begin(), result.end());
    return result;
}

const SimplexList &Network::get_facets()
{
    SimplexList result{};
    if (!facets_.has_value() && !simplex_tree_.has_value())
    {
        calc_facets_interactions();
    }
    else if (!facets_.has_value() && simplex_tree_.has_value())
    {
        calc_facets_simplex_tree();
    }
    return *facets_;
}

std::vector<uint32_t> Network::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    if (simplex_tree_.has_value())
    {
        return calc_degree_sequence_simplex_tree(simplex_dimension, neighbor_dimension);
    }
    else
    {
        return calc_degree_sequence_interactions(simplex_dimension, neighbor_dimension);
    }
}

std::vector<uint32_t> Network::calc_degree_sequence_interactions(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    const auto &facets{select_higher_dimensional_simplices(get_facets(), simplex_dimension)};
    const auto &simplex_list{get_skeleton_simplices(facets, simplex_dimension)};
    SimplexSet simplices{simplex_list.begin(), simplex_list.end()};
    std::vector<uint32_t> degree_sequence{};
    degree_sequence.reserve(simplices.size());

    std::mutex mutex;
    std::atomic<uint32_t> counter{0U};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](auto &&simplex)
        {
            if (++counter % 100000 == 0 && simplices.size() > 1000000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rC++: Calculating degree sequence ..."
                          << counter << " / " << simplices.size();
            }
            // Container of vertices with which the simplex forms a simplex of neighbor dimension

            SimplexSet combinations_of_remaining_vertices;

            // iterate over all facets
            for (const auto &facet : facets)
            {
                if (simplex.is_face(facet))
                {
            std::cout << __LINE__ << std::endl;
                    const auto difference{facet - simplex}; // vertices of facet that are not in the simplex

            std::cout << __LINE__ << std::endl;
                    const auto skeleton{difference.get_skeleton(neighbor_dimension - simplex_dimension - 1)};

            std::cout << __LINE__ << std::endl;
                    std::copy(skeleton.begin(), skeleton.end(), std::inserter(combinations_of_remaining_vertices, combinations_of_remaining_vertices.end()));
            std::cout << __LINE__ << std::endl;
                }
            }
            std::lock_guard<std::mutex> lock_guard(mutex);
            degree_sequence.push_back(combinations_of_remaining_vertices.size());
        });

    return degree_sequence;
}

std::vector<uint32_t> Network::calc_degree_sequence_simplex_tree(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(simplex_tree_.has_value());
    std::vector<uint32_t> result{};
    const auto &simplices(get_simplices());
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            if (simplex_tree_->dimension(simplex) == static_cast<int32_t>(simplex_dimension))
            {
                const auto degree{simplex_tree_->cofaces_simplex_range(simplex, neighbor_dimension - simplex_dimension).size()};
                std::lock_guard<std::mutex> lock{mutex};
                result.push_back(degree);
            }
        });

    std::sort(result.begin(), result.end());
    return result;
}

void Network::calc_facets_interactions()
{
    std::cout << __LINE__ << std::endl;
    sort_interactions(true);
    std::mutex mutex;
    facets_ = SimplexList{};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        interactions_->begin(),
        interactions_->end(),
        [&](auto &&interaction)
        {
            const auto first_it{interactions_->begin() + (&interaction - &(interactions_->at(0)))};
            if (++counter % 10000 == 0 && interactions_->size() > 100000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rC++: Extracting facets ... " << counter << " / " << interactions_->size();
            }
            auto first_is_face{false};
            for (auto second_it{first_it + 1}; second_it < interactions_->end(); ++second_it)
            {
                first_is_face = false;
                if (interaction.is_face(*second_it))
                {
                    first_is_face = true;
                    break;
                }
            }
            if (!first_is_face)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                facets_->push_back(interaction);
            }
        });
    std::cout << __LINE__ << std::endl;
}

void Network::sort_interactions(const bool ascending)
{
    sort_simplices(*interactions_, ascending);
}

void Network::calc_facets_simplex_tree()
{
    assert(simplex_tree_.has_value());
    const auto &simplex_handles{get_simplices()};
    std::mutex mutex;

    facets_ = SimplexList{};
    std::for_each(
        execution_policy,
        simplex_handles.begin(),
        simplex_handles.end(),
        [&](auto &&simplex_handle)
        {
            const auto cofaces{simplex_tree_->cofaces_simplex_range(simplex_handle, 1U)};
            if (cofaces.size() == 0U)
            {
                const auto facet{Simplex{get_vertices(simplex_handle)}};
                std::lock_guard<std::mutex> lock_guard{mutex};
                facets_->push_back(facet);
            }
        });
}

void Network::set_facets(const ISimplexList &facets)
{
    std::cout << __LINE__ << std::endl;
    facets_ = create_simplices(facets);
    std::cout << __LINE__ << std::endl;
}

const SimplexList &Network::get_interactions() const
{
    assert(interactions_ && "No interactions");
    return *interactions_;
}

void Network::set_interactions(const ISimplexList &interactions)
{
    interactions_ = create_simplices(interactions);
}

ISimplexList Network::get_skeleton(const ISimplexList &simplices_in) const
{
    std::cout << __LINE__ << std::endl;
    const auto &simplices{create_simplices(simplices_in)};
    const auto skeleton_simplices{get_skeleton_simplices(simplices, max_dimension_)};
    const auto result{create_raw_simplices(skeleton_simplices)};
    std::cout << __LINE__ << std::endl;
    return result;
}

ISimplexList Network::get_skeleton_interface(const Dimension max_dimension)
{
    std::cout << __LINE__ << std::endl;
    if (simplex_tree_.has_value())
    {
        return convert_to_raw_simplices<>(simplex_tree_->skeleton_simplex_range(max_dimension));
    }
    const auto skeleton{get_skeleton_simplices(*interactions_, max_dimension)};
    std::cout << __LINE__ << std::endl;
    return create_raw_simplices(skeleton);
}

void Network::create_simplex_tree_from_interactions()
{
    std::cout << __LINE__ << std::endl;
    assert(interactions_.has_value());
    reset_simplex_tree();
    add_vertices(vertices_);
    add_simplices(*interactions_);
    std::cout << __LINE__ << std::endl;
}

template <typename Iterator>
ISimplexList Network::convert_to_raw_simplices(const Iterator &simplex_range)
{
    std::cout << __LINE__ << std::endl;
    ISimplexList result{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplex_range.begin(),
        simplex_range.end(),
        [&](const auto &simplex)
        {
            VertexList vertices{};
            vertices.reserve(simplex_tree_->dimension(simplex) + 1U);
            for (auto vertex : simplex_tree_->simplex_vertex_range(simplex))
            {
                vertices.push_back(vertex);
            }
            std::sort(vertices.begin(), vertices.end());
            std::lock_guard<std::mutex> lock{mutex};
            result.push_back(vertices);
        });
    std::cout << __LINE__ << std::endl;
    return result;
}

void Network::set_simplices(const ISimplexList &simplices)
{
    std::cout << __LINE__ << std::endl;
    reset_simplex_tree();
    facets_ = std::nullopt;
    add_simplices_interface(simplices);
    std::cout << __LINE__ << std::endl;
}

void Network::keep_only_vertices(const VertexList &vertices)
{
    std::cout << __LINE__ << std::endl;
    vertices_ = vertices;
    filter_interactions(vertices);
    if (simplex_tree_.has_value())
    {
        create_simplex_tree_from_interactions();
    }
    std::cout << __LINE__ << std::endl;
}

void Network::filter_simplex_tree(const VertexList &vertices)
{
    std::cout << __LINE__ << std::endl;
    SimplexTree filtered_simplex_tree{};
    const auto &simplices{get_simplices()};

    std::mutex mutex{};
    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            auto keep_simplex{true};
            for (auto vertex : simplex_tree_->simplex_vertex_range(simplex))
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
                std::lock_guard<std::mutex> lock{mutex};
                filtered_simplex_tree.insert_simplex_and_subfaces(simplex_tree_->simplex_vertex_range(simplex));
            }
        });

    reset_simplex_tree();
    simplex_tree_ = filtered_simplex_tree;
    std::cout << __LINE__ << std::endl;
}

void Network::filter_interactions(const VertexList &vertices)
{
    std::cout << __LINE__ << std::endl;
    SimplexList filtered_interactions{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        interactions_->begin(),
        interactions_->end(),
        [&](const auto &interaction)
        {
            auto keep_interaction{true};
            for (auto vertex : interaction.vertices())
            {
                const auto vertex_in_interaction{std::find(vertices.begin(), vertices.end(), vertex) != vertices.end()};
                if (!vertex_in_interaction)
                {
                    keep_interaction = false;
                    break;
                }
            }
            if (keep_interaction)
            {
                std::lock_guard<std::mutex> lock{mutex};
                filtered_interactions.push_back(interaction);
            }
        });
    interactions_ = filtered_interactions;
    std::cout << __LINE__ << std::endl;
}

Dimension Network::get_max_dimension() const
{
    return max_dimension_;
}

void Network::set_max_dimension(const Dimension max_dimension)
{
    max_dimension_ = max_dimension;
}
