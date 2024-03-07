#include <atomic>
#include <execution>
#include <mutex>

#include "network.h"
#include "numpy_cpp_conversion.h"
#include "simplex.h"

Network::Network(const VertexList &vertices, const ISimplexList &interactions)
    : vertices_{vertices},
      interactions_{create_simplices(interactions)},
      facets_{std::nullopt},
      simplex_tree_{std::nullopt}
{
    sort_simplices(interactions_, true);
}

void Network::reset()
{
    reset_simplicial_complex();
    facets_ = std::nullopt;
}

void Network::add_simplices(const SimplexList &simplices, const Dimension dimension)
{
    initialize_simplicial_complex_if_needed();
    const auto representable_simplices{get_skeleton_simplices(simplices, dimension)};
    const auto total{representable_simplices.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        std::execution::seq,
        representable_simplices.begin(),
        representable_simplices.end(),
        [&](const auto &simplex)
        {
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::cout << "\rInsert simplices ... " << counter << " / " << total;
            }
            simplex_tree_->insert_simplex_and_subfaces(simplex.vertices());
        });

    if (total > 10000U)
    {
        std::cout << "\rInsert simplices ... " << total << " / " << total;
    }
}

void Network::add_vertices(const VertexList &vertices)
{
    SimplexList simplices{};
    simplices.reserve(vertices.size());

    std::transform(
        vertices.begin(),
        vertices.end(),
        std::back_inserter(simplices),
        [](const auto vertex)
        { return Simplex{VertexList{vertex}}; });

    add_simplices(simplices, 0);
}

const VertexList &Network::get_vertices_interface() const
{
    return vertices_;
}

void Network::set_vertices(const VertexList &vertices)
{
    vertices_ = vertices;
}

ISimplexList Network::get_interactions_interface() const
{
    return create_raw_simplices(get_interactions());
}

ISimplexList Network::get_facets_interface()
{
    return create_raw_simplices(get_facets());
}

ISimplexList Network::get_simplices_interface()
{
    return convert_to_raw_simplices(get_simplices());
}

template <typename Iterator>
ISimplexList Network::convert_to_raw_simplices(const Iterator &simplex_range)
{
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
    return result;
}

uint32_t Network::num_vertices()
{
    return vertices_.size();
}

std::vector<Dimension> Network::calc_interaction_dimension_distribution()
{
    const auto result{calc_dimension_distribution(interactions_)};
    return result;
}

std::vector<Dimension> Network::calc_facet_dimension_distribution()
{
    const auto &facets{get_facets()};
    const auto result{calc_dimension_distribution(facets)};
    return result;
}

const SimplexList &Network::get_facets()
{
    SimplexList result{};
    if (!facets_.has_value())
    {
        calc_facets();
    }
    return *facets_;
}

std::vector<uint32_t> Network::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    auto result{is_valid() ? calc_degree_sequence_simplicial_complex(simplex_dimension, neighbor_dimension)
                           : calc_degree_sequence_interactions(simplex_dimension, neighbor_dimension)};
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<uint32_t> Network::calc_degree_sequence_simplicial_complex(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(is_valid());
    assert(neighbor_dimension > simplex_dimension);
    std::vector<uint32_t> result{};
    // const auto &simplices(get_simplices());
    const auto &simplices(simplex_tree_->skeleton_simplex_range(simplex_dimension));
    std::mutex mutex{};
    // const auto total{simplices.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            if (++counter % 1000 == 0)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rCalc degree sequence (simplex tree) ... " << counter;
            }
            if (simplex_tree_->dimension(simplex) == simplex_dimension)
            {
                const auto degree{simplex_tree_->cofaces_simplex_range(simplex, neighbor_dimension - simplex_dimension).size()};
                std::lock_guard<std::mutex> lock{mutex};
                result.push_back(degree);
            }
        });

    std::cout << "\rCalc degree sequence (simplex tree) ... " << counter;

    return result;
}

std::vector<uint32_t> Network::calc_degree_sequence_interactions(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(neighbor_dimension > simplex_dimension);
    const auto &facets{select_higher_dimensional_simplices(get_facets(), simplex_dimension)};
    const auto &simplex_list{get_skeleton_simplices(facets, simplex_dimension)};
    SimplexSet simplices{simplex_list.begin(), simplex_list.end()};
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
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rCalc degree sequence (interactions) ... "
                          << counter << " / " << total;
            }
            // Container of vertices with which the simplex forms a simplex of neighbor dimension

            SimplexSet combinations_of_remaining_vertices;

            // iterate over all facets
            std::for_each(
                execution_policy,
                facets.begin(),
                facets.end(),
                [&](const auto &facet)
                {
                    if (simplex.is_face(facet))
                    {
                        const auto difference{facet - simplex}; // vertices of facet that are not in the simplex
                        const auto skeleton{difference.get_skeleton(neighbor_dimension - simplex_dimension - 1)};

                        std::lock_guard<std::mutex> lock_guard(mutex);
                        std::copy(skeleton.begin(), skeleton.end(), std::inserter(combinations_of_remaining_vertices, combinations_of_remaining_vertices.end()));
                    }
                });
            const auto degree{combinations_of_remaining_vertices.size()};
            std::lock_guard<std::mutex> lock_guard(mutex);
            degree_sequence.push_back(degree);
        });

    if (total > 10000U)
    {
        std::cout << "\rCalc degree sequence (interactions) ... " << total << " / " << total;
    }
    return degree_sequence;
}

void Network::calc_facets()
{
    std::mutex mutex{};
    facets_ = SimplexList{};
    const auto total{interactions_.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        interactions_.begin(),
        interactions_.end(),
        [&](auto &&interaction)
        {
            const auto first_it{interactions_.begin() + (&interaction - &(interactions_[0]))};
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rCalc facets (interactions) ... " << counter << " / " << total;
            }
            auto first_is_face{false};
            for (auto second_it{first_it + 1}; second_it < interactions_.end(); ++second_it)
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

    if (total > 10000U)
    {
        std::cout << "\rCalc facets (interactions) ... " << total << " / " << total;
    }
}

const SimplexList &Network::get_interactions() const
{
    return interactions_;
}

void Network::set_interactions(const ISimplexList &interactions)
{
    interactions_ = create_simplices(interactions);
}

ISimplexList Network::get_skeleton_interface(const Dimension max_dimension)
{
    return create_raw_simplices(get_skeleton(max_dimension));
}

void Network::create_simplicial_complex(const Dimension max_dimension)
{
    reset_simplicial_complex();
    add_vertices(vertices_);
    add_simplices(interactions_, max_dimension);
    facets_ = std::nullopt;
}

void Network::keep_only_vertices(const VertexList &vertices)
{
    vertices_ = vertices;
    filter_interactions(vertices);
    if (is_valid())
    {
        const auto max_dimension{simplex_tree_->dimension()};
        create_simplicial_complex(max_dimension);
    }
}

SimplexList Network::get_skeleton(const Dimension max_dimension)
{
    if (is_valid())
    {
        return get_skeleton_simplicial_complex(max_dimension);
    }
    return get_skeleton_interactions(max_dimension);
}

void Network::filter_interactions(const VertexList &vertices)
{
    SimplexList filtered_interactions{};
    std::mutex mutex{};
    std::for_each(
        execution_policy,
        interactions_.begin(),
        interactions_.end(),
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
}

void Network::expand(const Dimension max_dimension)
{
    assert(is_valid());
    simplex_tree_->expansion(max_dimension + 1U);
}

void Network::initialize_simplicial_complex_if_needed()
{
    if (!is_valid())
    {
        simplex_tree_ = SimplexTree{};
    }
}

bool Network::is_valid() const
{
    return simplex_tree_.has_value();
}

void Network::reset_simplicial_complex()
{
    simplex_tree_ = std::nullopt;
}

std::vector<Dimension> Network::calc_simplex_dimension_distribution()
{
    assert(is_valid());
    const auto &simplices{get_simplices()};
    std::vector<Dimension> result{};

    std::mutex mutex{};
    const auto total{simplices.size()};
    std::atomic<uint32_t> counter{0U};
    result.reserve(total);

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            const auto dimension{simplex_tree_->dimension(simplex)};

            std::lock_guard<std::mutex> lock{mutex};
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::cout << "\rCalc dimension distribution (simplex tree) ... "
                          << counter << " / " << total;
            }
            result.push_back(dimension);
        });

    if (total > 10000U)
    {
        std::cout << "\rCalc dimension distribution (simplex tree) ... " << total << " / " << total;
    }

    std::sort(result.begin(), result.end());
    return result;
}

VertexList Network::get_vertices(const SimplexHandle &simplex_handle)
{
    assert(is_valid());
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