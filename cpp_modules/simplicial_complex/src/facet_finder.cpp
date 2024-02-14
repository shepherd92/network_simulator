#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#include "combination.h"
#include "facet_finder.h"
#include "simplex.h"

void combinations(
    const std::vector<int32_t> &elements,
    const uint32_t k,
    std::set<std::vector<int32_t>> &subarrays,
    std::vector<int32_t> &out,
    const uint32_t i);

std::vector<std::vector<int32_t>> extract_facets(const std::vector<std::vector<int32_t>> &simplices_in)
{
  auto simplices{sort_simplices(create_simplices(simplices_in), true)};
  std::vector<std::vector<int32_t>> facets;
  auto counter{0U};
  for (auto first{simplices.begin()}; first < simplices.end(); ++first)
  {
    if (++counter % 10000 == 0 && simplices.size() > 100000U)
    {
      std::cout << "\rC++: Extracting facets ... " << counter << " / " << simplices.size();
    }

    auto first_is_face{false};
    for (auto second{first + 1}; second < simplices.end(); ++second)
    {
      first_is_face = false;
      if (first->is_face(*second))
      {
        first_is_face = true;
        break;
      }
    }
    if (!first_is_face)
    {
      std::vector<int32_t> vertices_vector;
      vertices_vector.assign(first->vertices().begin(), first->vertices().end());
      facets.push_back(vertices_vector);
    }
  }

  return facets;
}
