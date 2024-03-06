#include <atomic>
#include <execution>
#include <mutex>

#include "network.h"
#include "numpy_cpp_conversion.h"
#include "simplex.h"

Network::Network()
    : vertices_{},
      interactions_{std::nullopt},
      facets_{std::nullopt}
{
}

void Network::reset()
{
    reset_simplicial_complex();
    interactions_ = std::nullopt;
    facets_ = std::nullopt;
}

void Network::add_simplices_interface(const ISimplexList &simplices_in)
{
    initialize_simplicial_complex_if_needed();
    const auto simplices{create_simplices(simplices_in)};
    add_simplices(simplices);
    reset_persistence();
}

void Network::add_simplices(const SimplexList &simplices)
{
    initialize_simplicial_complex_if_needed();
    const auto representable_simplices{convert_to_representable_simplices(simplices)};
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
            add_simplex(simplex.vertices());
        });

    if (total > 10000U)
    {
        std::cout << "\rInsert simplices ... " << total << " / " << total;
    }
}

void Network::add_vertices(const VertexList &vertices)
{
    ISimplexList simplices{};
    simplices.reserve(vertices.size());

    std::transform(
        vertices.begin(),
        vertices.end(),
        std::back_inserter(simplices),
        [](const auto vertex)
        { return VertexList{vertex}; });

    add_simplices_interface(simplices);

    reset_persistence();
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

uint32_t Network::num_vertices()
{
    return vertices_.size();
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

const SimplexList &Network::get_facets()
{
    SimplexList result{};
    if (!facets_.has_value() && !is_valid())
    {
        calc_facets_interactions();
    }
    else if (!facets_.has_value() && is_valid())
    {
        calc_facets_simplicial_complex();
    }
    assert(facets_.has_value() && "Facets not calculated");
    return *facets_;
}

std::vector<uint32_t> Network::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    if (is_valid())
    {
        return calc_degree_sequence_simplicial_complex(simplex_dimension, neighbor_dimension);
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

void Network::calc_facets_interactions()
{
    sort_interactions(true);
    std::mutex mutex{};
    facets_ = SimplexList{};
    const auto total{interactions_->size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        interactions_->begin(),
        interactions_->end(),
        [&](auto &&interaction)
        {
            const auto first_it{interactions_->begin() + (&interaction - &(interactions_->at(0)))};
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rCalc facets (interactions) ... " << counter << " / " << total;
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

    if (total > 10000U)
    {
        std::cout << "\rCalc facets (interactions) ... " << total << " / " << total;
    }
}

void Network::sort_interactions(const bool ascending)
{
    sort_simplices(*interactions_, ascending);
}

void Network::set_facets(const ISimplexList &facets)
{
    facets_ = create_simplices(facets);
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

ISimplexList Network::get_skeleton_interface(const Dimension max_dimension)
{
    return create_raw_simplices(get_skeleton(max_dimension));
}

SimplexList Network::get_skeleton(const Dimension max_dimension)
{
    if (is_valid())
    {
        return get_skeleton_simplicial_complex(max_dimension);
    }
    return get_skeleton_simplices(*interactions_, max_dimension);
}

void Network::create_simplicial_complex_from_interactions()
{
    assert(interactions_.has_value());
    reset_simplicial_complex();
    add_vertices(vertices_);
    add_simplices(*interactions_);
}

template <typename Iterator>
ISimplexList Network::convert_to_raw_simplices(const Iterator &simplex_range)
{
    ISimplexList result{};
    std::mutex mutex{};
    std::atomic<uint32_t> counter{0U};

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
            if (++counter % 1000 == 0)
            {
                std::cout << "\rConvert to raw simplices ... " << counter;
            }
            result.push_back(vertices);
        });

    return result;
}

void Network::keep_only_vertices(const VertexList &vertices)
{
    vertices_ = vertices;
    filter_interactions(vertices);
    if (is_valid())
    {
        create_simplicial_complex_from_interactions();
    }
}

void Network::filter_interactions(const VertexList &vertices)
{
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
}
