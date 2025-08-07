#include "finite_network.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

FiniteNetwork::FiniteNetwork()
    : simplex_tree_{std::nullopt},
      persistent_cohomology_{nullptr}
{
}

FiniteNetwork::~FiniteNetwork()
{
    reset_persistence();
}

SimplexList FiniteNetwork::get_skeleton(const Dimension max_dimension)
{
    SimplexList result{};
    for (auto dimension{0}; dimension <= max_dimension; ++dimension)
    {
        result += get_simplices(dimension);
    }
    return result;
}

void FiniteNetwork::reset()
{
    Network::reset();
    reset_simplicial_complex();
}

void FiniteNetwork::reset_simplicial_complex()
{
    simplex_tree_ = std::nullopt;
    reset_persistence();
}

void FiniteNetwork::reset_persistence()
{
    delete persistent_cohomology_;
    persistent_cohomology_ = nullptr;
}

FiniteNetwork::PersistentCohomology &FiniteNetwork::get_persistence()
{
    if (!persistent_cohomology_)
    {
        calc_persistent_cohomology();
    }
    return *persistent_cohomology_;
}

void FiniteNetwork::calc_persistent_cohomology()
{
    reset_persistence();
    assert_simplicial_complex_is_built();
    std::cout << "\rCompute persistent cohomology..." << std::flush;
    persistent_cohomology_ = new PersistentCohomology{*simplex_tree_};
    persistent_cohomology_->init_coefficients(2);
    persistent_cohomology_->compute_persistent_cohomology();
    std::cout << "done" << std::flush;
}

std::vector<uint32_t> FiniteNetwork::calc_coface_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(neighbor_dimension > simplex_dimension);

    const auto &cofaces{get_simplices(neighbor_dimension)};
    auto simplex_degree_map{cofaces.calc_degree_sequence(simplex_dimension)};

    // order of the degree values does not matter
    std::vector<uint32_t> result{};
    const auto &simplices{get_simplices(simplex_dimension)};
    result.reserve(simplex_degree_map.size());
    for (const auto &simplex : simplices)
    {
        result.emplace_back(simplex_degree_map[simplex]);
    }

    return result;
}

std::vector<int32_t> FiniteNetwork::calc_betti_numbers()
{
    std::vector<int32_t> result{max_dimension_, 0};
    if (num_simplices(0) == 1U)
    {
        // handle an error in Gudhi
        result[0] = 1;
    }
    else
    {
        const auto &persistent_cohomology{get_persistence()};
        result = persistent_cohomology.betti_numbers();
        result.resize(max_dimension_);
    }
    assert(static_cast<int32_t>(result.size()) == max_dimension_);
    return result;
}

void FiniteNetwork::assert_simplicial_complex_is_built()
{
    if (!is_valid())
    {
        create_simplicial_complex();
    }
}

void FiniteNetwork::assert_simplicial_complex_is_initialized()
{
    if (!is_valid())
    {
        simplex_tree_ = SimplexTree{};
    }
}

bool FiniteNetwork::is_valid() const
{
    return simplex_tree_.has_value();
}

void FiniteNetwork::create_simplicial_complex()
{
    reset_simplicial_complex();
    add_vertices();
    fill_simplicial_complex();
}

void FiniteNetwork::add_vertices()
{
    assert_simplicial_complex_is_initialized();
    get_simplex_tree()->insert_batch_vertices(get_vertices());
}

PointIdList FiniteNetwork::get_simplex_vertices(const SimplexHandle &simplex_handle)
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

std::optional<FiniteNetwork::SimplexTree> &FiniteNetwork::get_simplex_tree()
{
    return simplex_tree_;
}