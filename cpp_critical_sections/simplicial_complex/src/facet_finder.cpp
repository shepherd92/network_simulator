#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "facet_finder.h"

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
  auto counter{0U};

  for (auto first = sorted_simplices.begin(); first != sorted_simplices.end(); ++first)
  {
    if (++counter % 1000 == 0)
    {
      std::cout << "\rC++: extracting facets..."
                << std::setprecision(3)
                << static_cast<float>(counter) / static_cast<float>(sorted_simplices.size()) * 100
                << "%    ";
    }

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

std::vector<std::vector<int32_t>> select_simplices_by_dimension(
    const std::vector<std::vector<int32_t>> &simplices,
    const uint32_t dimension)
{
  std::vector<std::vector<int32_t>> selected_simplices;

  for (const auto &simplex : simplices)
  {
    if (simplex.size() == dimension + 1)
    {
      selected_simplices.push_back(simplex);
    }
  }

  return sort_simplices(selected_simplices);
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

  std::vector<int32_t> degree_sequence;
  auto counter{0U};

  for (const auto &simplex : selected_simplices)
  {
    if (++counter % 1000 == 0)
    {
      std::cout << "\rC++: calculating degree sequence..."
                << std::setprecision(3)
                << static_cast<float>(counter) / static_cast<float>(selected_simplices.size()) * 100.
                << "%    ";
    }

    int higher_order_degree{0};
    for (const auto &facet : possible_neighbors)
    {
      if (includes(facet.begin(), facet.end(), simplex.begin(), simplex.end()))
      {
        higher_order_degree += choose(facet.size() - 1 - simplex_dimension, neighbor_dimension - simplex_dimension);
      }
    }
    degree_sequence.push_back(higher_order_degree);
  }

  return degree_sequence;
}
