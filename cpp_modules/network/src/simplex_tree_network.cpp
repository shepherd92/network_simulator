#include "simplex.h"
#include "simplex_tree_network.h"

SimplexTreeNetwork::SimplexTreeNetwork(const Dimension max_dimension)
    : Network{},
      max_dimension_{max_dimension},
      simplex_tree_{std::nullopt}
{
}

void SimplexTreeNetwork::initialize_simplicial_complex_if_needed()
{
    if (!is_valid())
    {
        simplex_tree_ = SimplexTree{};
    }
}

bool SimplexTreeNetwork::is_valid() const
{
    return simplex_tree_.has_value();
}

void SimplexTreeNetwork::reset_simplicial_complex()
{
    simplex_tree_ = std::nullopt;
}

void SimplexTreeNetwork::expand()
{
    assert(simplex_tree_.has_value());
    simplex_tree_->expansion(max_dimension_ + 1U);
}

void SimplexTreeNetwork::add_simplex(const VertexList &simplex)
{
    simplex_tree_->insert_simplex_and_subfaces(simplex);
}

SimplexTreeNetwork::SimplexTree &SimplexTreeNetwork::get_simplex_tree()
{
    assert(is_valid() || interactions_.has_value());
    if (!is_valid())
    {
        create_simplicial_complex_from_interactions();
    }
    return *simplex_tree_;
}

std::vector<uint32_t> SimplexTreeNetwork::calc_degree_sequence_simplicial_complex(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(is_valid());
    assert(neighbor_dimension > simplex_dimension);

    std::vector<uint32_t> result{};
    const auto &simplices(get_simplices());
    std::mutex mutex{};
    const auto total{simplices.size()};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const auto &simplex)
        {
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rCalc degree sequence (simplex tree) ... " << counter << " / " << total;
            }
            if (simplex_tree_->dimension(simplex) == simplex_dimension)
            {
                const auto degree{simplex_tree_->cofaces_simplex_range(simplex, neighbor_dimension - simplex_dimension).size()};
                std::lock_guard<std::mutex> lock{mutex};
                result.push_back(degree);
            }
        });

    if (total > 10000U)
    {
        std::cout << "\rCalc degree sequence (simplex tree) ... " << total << " / " << total;
    }

    std::sort(result.begin(), result.end());
    return result;
}

std::vector<Dimension> SimplexTreeNetwork::calc_simplex_dimension_distribution()
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

SimplexList SimplexTreeNetwork::convert_to_representable_simplices(const SimplexList &simplices_in) const
{
    const auto simplices{get_skeleton_simplices(simplices_in, max_dimension_)};
    return simplices;
}

void SimplexTreeNetwork::calc_facets_simplicial_complex()
{
    assert(is_valid());
    const auto &simplex_handles{get_simplices()};
    std::mutex mutex{};
    const auto total{simplex_handles.size()};
    std::atomic<uint32_t> counter{0U};

    facets_ = SimplexList{};
    std::for_each(
        execution_policy,
        simplex_handles.begin(),
        simplex_handles.end(),
        [&](auto &&simplex_handle)
        {
            if (++counter % 1000 == 0 && total > 10000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                std::cout << "\rCalc facets (simplex tree) ... " << counter << " / " << total;
            }
            const auto cofaces{simplex_tree_->cofaces_simplex_range(simplex_handle, 1U)};
            if (cofaces.size() == 0U)
            {
                std::lock_guard<std::mutex> lock_guard{mutex};
                facets_->emplace_back(Simplex{get_vertices(simplex_handle)});
            }
        });

    if (total > 10000U)
    {
        std::cout << "\rCalc facets (simplex tree) ... " << total << " / " << total;
    }
}

SimplexList SimplexTreeNetwork::get_skeleton_simplicial_complex(const Dimension max_dimension)
{
    assert(is_valid());
    const auto &simplices{simplex_tree_->skeleton_simplex_range(max_dimension)};
    SimplexList skeleton_simplices{};
    std::transform(
        simplices.begin(),
        simplices.end(),
        std::back_inserter(skeleton_simplices),
        [this](const auto &simplex_handle)
        {
            return Simplex{get_vertices(simplex_handle)};
        });
    return skeleton_simplices;
}

VertexList SimplexTreeNetwork::get_vertices(const SimplexHandle &simplex_handle)
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

Dimension SimplexTreeNetwork::get_max_dimension() const
{
    return max_dimension_;
}

void SimplexTreeNetwork::set_max_dimension(const Dimension max_dimension)
{
    max_dimension_ = max_dimension;
}
