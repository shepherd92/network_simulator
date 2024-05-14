#include <atomic>
#include <cassert>
#include <execution>
#include <mutex>

#include "network.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"
#include "typedefs.h"

Network::Network(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const SimplexList &nonempty_interactions,
    const uint32_t num_of_empty_interactions)
    : max_dimension_{max_dimension},
      vertices_{vertices},
      nonempty_interactions_{nonempty_interactions},
      num_of_empty_interactions_{num_of_empty_interactions},
      facets_{std::nullopt},
      simplices_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
    if (vertices_.empty())
    {
        vertices_ = nonempty_interactions_.vertices();
    }
    std::sort(vertices_.begin(), vertices_.end());
}

void Network::reset()
{
    facets_ = std::nullopt;
    simplices_ = std::vector<std::optional<SimplexList>>{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt};
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
    return get_interactions().raw();
}

ISimplexList Network::get_facets_interface()
{
    return get_facets().raw();
}

uint32_t Network::num_vertices()
{
    return vertices_.size();
}

std::vector<Dimension> Network::calc_interaction_dimension_distribution() const
{
    return nonempty_interactions_.calc_dimension_distribution();
}

std::vector<Dimension> Network::calc_facet_dimension_distribution()
{
    return get_facets().calc_dimension_distribution();
}

const SimplexList &Network::get_facets()
{
    if (!facets_.has_value())
    {
        facets_ = nonempty_interactions_.facets();
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
    const auto &possible_cofaces{get_neighbors(neighbor_dimension)};
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

const SimplexList &Network::get_interactions() const
{
    return nonempty_interactions_;
}

void Network::set_interactions(const ISimplexList &interactions)
{
    nonempty_interactions_ = SimplexList{interactions};
}

ISimplexList Network::get_skeleton_interface(const Dimension max_dimension)
{
    return get_skeleton(max_dimension).raw();
}

void Network::keep_only_vertices(const PointIdList &vertices)
{
    vertices_ = vertices;
    nonempty_interactions_ = nonempty_interactions_.filter(vertices);
    reset();
}

SimplexList Network::get_skeleton(const Dimension max_dimension)
{
    SimplexList result{};
    for (auto dimension{0}; dimension <= max_dimension; ++dimension)
    {
        result += get_simplices(dimension);
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

uint32_t Network::num_simplices(const Dimension dimension)
{
    return get_simplices(dimension).size();
}