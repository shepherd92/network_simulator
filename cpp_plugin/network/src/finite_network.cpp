#include <algorithm>

#include "finite_network.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"

FiniteNetwork::FiniteNetwork()
{
}

FiniteNetwork::FiniteNetwork(FiniteNetwork &&) noexcept
{
}

FiniteNetwork &FiniteNetwork::operator=(FiniteNetwork &&) noexcept
{
    return *this;
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
    std::vector<int32_t> result(max_dimension_, 0);
    if (num_simplices(0) == 1U)
    {
        // handle an error in Gudhi
        result[0] = 1;
        return result;
    }
    else
    {
        assert_persistence_cohomology_is_calculated();
        const auto betti_numbers{persistent_cohomology_->betti_numbers()};
        for (auto dimension{0};
             dimension < std::min(max_dimension_, static_cast<int32_t>(betti_numbers.size()));
             ++dimension)
        {
            result[dimension] = betti_numbers[dimension];
        }
    }
    return result;
}
