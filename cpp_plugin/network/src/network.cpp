#include <atomic>
#include <cassert>
#include <execution>
#include <mutex>

#include "network.h"
#include "simplex_list.h"
#include "tools.h"
#include "typedefs.h"

Network::Network(
    const Dimension max_dimension,
    const PointIdList &vertices)
    : max_dimension_{max_dimension},
      simplices_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt},
      vertices_{vertices}
{
    std::sort(vertices_.begin(), vertices_.end());
}

void Network::reset()
{
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

PointIdList Network::get_vertices() const
{
    return vertices_;
}

uint32_t Network::num_vertices()
{
    return vertices_.size();
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

ISimplexList Network::get_skeleton_interface(const Dimension max_dimension)
{
    return get_skeleton(max_dimension).raw();
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

void Network::keep_only_vertices(const PointIdList &vertices)
{
    set_vertices(vertices);
    reset();
}