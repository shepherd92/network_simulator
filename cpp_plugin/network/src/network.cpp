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
      simplex_tree_{std::nullopt},
      persistent_cohomology_{nullptr},
      simplices_{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt}
{
    std::sort(vertices_.begin(), vertices_.end());
}

Network::~Network()
{
    reset_persistence();
}

Network::Network(Network &&other) noexcept
{
    max_dimension_ = std::move(other.max_dimension_);
    vertices_ = std::move(other.vertices_);
    simplex_tree_ = std::move(other.simplex_tree_);
    persistent_cohomology_ = other.persistent_cohomology_;
    simplices_ = std::move(other.simplices_);
}

Network &Network::operator=(Network &&other) noexcept
{
    if (this != &other)
    {
        max_dimension_ = std::move(other.max_dimension_);
        vertices_ = std::move(other.vertices_);
        simplex_tree_ = std::move(other.simplex_tree_);
        persistent_cohomology_ = other.persistent_cohomology_;
        simplices_ = std::move(other.simplices_);
    }
    return *this;
}

void Network::reset()
{
    reset_simplicial_complex();
    simplices_ = std::vector<std::optional<SimplexList>>{static_cast<uint32_t>(max_dimension_) + 1U, std::nullopt};
}

void Network::reset_simplicial_complex()
{
    simplex_tree_ = std::nullopt;
    reset_persistence();
}

void Network::reset_persistence()
{
    delete persistent_cohomology_;
    persistent_cohomology_ = nullptr;
}

Network::PersistentCohomology &Network::get_persistence()
{
    if (!persistent_cohomology_)
    {
        calc_persistent_cohomology();
    }
    return *persistent_cohomology_;
}

void Network::calc_persistent_cohomology()
{
    reset_persistence();
    assert_simplicial_complex_is_built();
    std::cout << "\rCompute persistent cohomology..." << std::flush;
    persistent_cohomology_ = new PersistentCohomology{*simplex_tree_};
    persistent_cohomology_->init_coefficients(2);
    persistent_cohomology_->compute_persistent_cohomology();
    std::cout << "done" << std::flush;
}

void Network::assert_simplicial_complex_is_built()
{
    if (!is_valid())
    {
        create_simplicial_complex();
    }
}

void Network::assert_simplicial_complex_is_initialized()
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

void Network::create_simplicial_complex()
{
    add_vertices();
    fill_simplicial_complex();
}

void Network::add_vertices()
{
    assert_simplicial_complex_is_initialized();
    simplex_tree_->insert_batch_vertices(vertices_);
}

PointIdList Network::get_simplex_vertices(const SimplexHandle &simplex_handle)
{
    assert_simplicial_complex_is_built();
    PointIdList result{};
    if (simplex_handle != simplex_tree_->null_simplex())
    {
        for (const auto &vertex : simplex_tree_->simplex_vertex_range(simplex_handle))
        {
            result.push_back(vertex);
        }
    }
    return result;
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
    std::vector<Dimension> result;
    for (auto dimension{0}; dimension <= max_dimension_; ++dimension)
    {
        std::vector<Dimension> result_for_dimension(get_simplices(dimension).size(), dimension);
        result.insert(result.end(), result_for_dimension.begin(), result_for_dimension.end());
    }
    return result;
}
