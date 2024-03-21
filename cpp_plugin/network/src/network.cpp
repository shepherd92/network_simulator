#include <atomic>
#include <execution>
#include <mutex>

#include "network.h"
#include "simplex.h"
#include "tools.h"
#include "typedefs.h"

Network::Network(const Dimension max_dimension, const PointIdList &vertices, const SimplexList &interactions)
    : max_dimension_{max_dimension},
      vertices_{vertices},
      interactions_{interactions},
      facets_{std::nullopt},
      simplices_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
    sort_simplices(interactions_, true);
    if (vertices_.empty())
    {
        std::unordered_set<PointId> unique_vertices{};
        std::for_each(
            std::execution::seq,
            interactions_.begin(),
            interactions_.end(),
            [&](const auto &interaction)
            {
                for (const auto &vertex : interaction.vertices())
                {
                    unique_vertices.insert(vertex);
                }
            });
        vertices_ = PointIdList{unique_vertices.begin(), unique_vertices.end()};
    }
    std::sort(vertices_.begin(), vertices_.end());
}

void Network::reset()
{
    facets_ = std::nullopt;
}

Dimension Network::get_max_dimension() const
{
    return max_dimension_;
}

void Network::set_max_dimension(const Dimension dimension)
{
    max_dimension_ = dimension;
}

void Network::set_vertices(const PointIdList &vertices)
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

std::vector<Dimension> Network::calc_interaction_dimension_distribution() const
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
    if (!facets_.has_value())
    {
        facets_ = calc_facets();
    }
    return *facets_;
}

const SimplexList &Network::get_simplices(const Dimension dimension)
{
    assert(dimension <= max_dimension_);
    if (!simplices_[dimension].has_value())
    {
        simplices_[dimension] = calc_simplices(dimension);
    }
    return *simplices_[dimension];
}

std::vector<uint32_t> Network::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(neighbor_dimension > simplex_dimension);
    const auto &simplices{get_simplices(simplex_dimension)};
    const auto &possible_cofaces{get_simplices(neighbor_dimension)};
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

SimplexList Network::calc_facets() const
{
    std::mutex mutex{};
    SimplexList facets{};
    const auto total{interactions_.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        interactions_.begin(),
        interactions_.end(),
        [&](auto &&interaction)
        {
            const auto first_it{interactions_.begin() + (&interaction - &(interactions_[0]))};
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
                facets.push_back(interaction);
            }
            log_progress(++counter, total, 1000U, "Calc facets");
        });

    log_progress(counter, total, 1U, "Calc facets");
    return facets;
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

void Network::keep_only_vertices(const PointIdList &vertices)
{
    vertices_ = vertices;
    interactions_ = filter_simplices(interactions_, vertices);
    reset();
}

SimplexList Network::get_skeleton(const Dimension max_dimension)
{
    SimplexList result{};
    for (auto dimension{0}; dimension <= max_dimension; ++dimension)
    {
        const auto simplices_of_dimension{get_simplices(dimension)};
        result.reserve(result.size() + simplices_of_dimension.size());
        result.insert(result.end(), simplices_of_dimension.begin(), simplices_of_dimension.end());
    }
    return result;
}

std::vector<Dimension> Network::calc_simplex_dimension_distribution()
{
    std::vector<Dimension> result{max_dimension_ + 1, 0};
    for (auto dimension{0}; dimension <= max_dimension_; ++dimension)
    {
        result[dimension] = get_simplices(dimension).size();
    }
    return result;
}

std::vector<uint32_t> Network::calc_vertex_interaction_degree_distribution() const
{
    const auto vertices{get_vertices()};
    std::vector<uint32_t> result{};
    result.reserve(vertices.size());
    std::for_each(
        std::execution::seq,
        vertices.begin(),
        vertices.end(),
        [&](const auto &vertex)
        {
            result.emplace_back(get_cofaces(interactions_, Simplex{PointIdList{vertex}}).size());
        });
    return result;
}

uint32_t Network::num_simplices(const Dimension dimension)
{
    return get_simplices(dimension).size();
}