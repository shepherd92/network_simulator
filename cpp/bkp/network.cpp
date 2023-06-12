#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#include "network.h"

void combinations(
    const std::vector<int32_t> &elements,
    const uint32_t k,
    std::set<std::vector<int32_t>> &subarrays,
    std::vector<int32_t> &out,
    const uint32_t i);

std::vector<simplex> Network::get_facets()
{
  return facets;
}

std::vector<std::vector<int32_t>> sort_simplices(const std::vector<std::vector<int32_t>> &simplices)
{
  auto unique_simplices{simplices};
  unique_simplices.erase(unique(unique_simplices.begin(), unique_simplices.end()), unique_simplices.end());

  // sort list of simplices
  std::sort(unique_simplices.begin(), unique_simplices.end(), [](const std::vector<int32_t> &a, const std::vector<int32_t> &b)
            { return a.size() < b.size(); });

  // sort simplices themselves
  for (auto &simplex : unique_simplices)
  {
    std::sort(simplex.begin(), simplex.end());
  }

  return unique_simplices;
}

uint32_t choose(const uint32_t n, uint32_t k)
{
  if (k > n)
    return 0;
  if (k * 2 > n)
    k = n - k;
  if (k == 0)
    return 1;

  auto result{n};
  for (auto i{2U}; i <= k; ++i)
  {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

std::vector<std::vector<int32_t>> extract_facets(const std::vector<std::vector<int32_t>> &simplices)
{
  auto sorted_simplices{sort_simplices(simplices)};
  std::vector<std::vector<int32_t>> facets;

  for (auto first = sorted_simplices.begin(); first != sorted_simplices.end(); ++first)
  {
    auto first_is_subset{false};
    for (auto second{next(first)}; second != sorted_simplices.end(); ++second)
    {
      first_is_subset = false;
      if (includes(second->begin(), second->end(), first->begin(), first->end()))
      {
        first_is_subset = true;
        break;
      }
    }
    if (!first_is_subset)
    {
      facets.push_back(*first);
    }
  }

  return facets;
}

std::vector<int32_t> calc_degree_sequence(
    const std::vector<std::vector<int32_t>> &simplices,
    const std::vector<std::vector<int32_t>> &facets,
    const uint32_t simplex_dimension,
    const uint32_t neighbor_dimension)
{
  const auto selected_simplices{select_simplices_by_dimension(simplices, simplex_dimension)};

  std::vector<std::vector<int32_t>> possible_neighbors;
  copy_if(facets.begin(), facets.end(),
          back_inserter(possible_neighbors),
          [neighbor_dimension](const std::vector<int32_t> &facet)
          { return facet.size() > neighbor_dimension; });

  const auto sorted_possible_neighbors{sort_simplices(possible_neighbors)};

  std::vector<int32_t> degree_sequence;

  for (const auto &simplex : selected_simplices)
  {
    // Container of vertices with which the simplex forms a simplex of neighbor dimension
    std::set<std::vector<int32_t>> combinations_of_remaining_vertices;

    // iterate over all facets
    for (const auto &facet : sorted_possible_neighbors)
    {
      // check if the facet includes the vertices of the simplex
      if (includes(facet.begin(), facet.end(), simplex.begin(), simplex.end()))
      {
        std::vector<int32_t> remaining_vertices; // vertices of facet that are not in the simplex
        std::set_difference(
            facet.begin(), facet.end(),
            simplex.begin(), simplex.end(),
            std::back_inserter(remaining_vertices));
        std::sort(remaining_vertices.begin(), remaining_vertices.end());

        std::vector<int32_t> temp_combination;

        combinations(
            remaining_vertices,
            neighbor_dimension - simplex_dimension,
            combinations_of_remaining_vertices,
            temp_combination,
            0U);
      }
    }
    degree_sequence.push_back(combinations_of_remaining_vertices.size());
  }

  return degree_sequence;
}