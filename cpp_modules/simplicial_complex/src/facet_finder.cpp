#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#include "facet_finder.h"
#include "simplex.h"

std::vector<VertexList> extract_facets_interface(const std::vector<VertexList> &simplices_in)
{
  auto simplices{sort_simplices(create_simplices(simplices_in), true)};
  std::vector<VertexList> facets;
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
      VertexList raw_simplex;
      raw_simplex.assign(first->vertices().begin(), first->vertices().end());
      facets.push_back(raw_simplex);
    }
  }

  return facets;
}
