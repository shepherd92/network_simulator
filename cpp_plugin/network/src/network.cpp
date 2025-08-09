#include <atomic>
#include <cassert>
#include <execution>
#include <mutex>

#include "network.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"
#include "typedefs.hpp"

Network::Network(
    const Dimension max_dimension,
    const PointIdList &vertices)
    : max_dimension_{max_dimension},
      vertices_{vertices},
      simplices_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
    std::sort(vertices_.begin(), vertices_.end());
}

void Network::reset()
{
    simplices_ = std::vector<std::optional<SimplexList>>{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt};
}

uint32_t Network::num_vertices()
{
    return vertices_.size();
}

const PointIdList &Network::get_vertices() const
{
    return vertices_;
}

void Network::set_vertices(const PointIdList &vertices)
{
    vertices_ = vertices;
    reset();
}

Dimension Network::get_max_dimension() const
{
    return max_dimension_;
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

uint32_t Network::num_simplices(const Dimension dimension)
{
    return get_simplices(dimension).size();
}

ISimplexList Network::get_skeleton_interface(const Dimension max_dimension)
{
    assert(max_dimension <= max_dimension_);
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
