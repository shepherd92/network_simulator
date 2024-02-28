#include <algorithm>
#include <execution>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <set>
#include <vector>

#include "facet_finder.h"
#include "simplex.h"

std::vector<VertexList> extract_facets_interface(const std::vector<VertexList> &simplices_in)
{
  auto simplices{sort_simplices(create_simplices(simplices_in), true)};
  std::mutex mutex;
  std::vector<VertexList> facets{};
  std::atomic<uint32_t> counter{0U};

  std::for_each(
      simplices.begin(),
      simplices.end(),
      [&](auto &&simplex)
      {
        const auto first_it{simplices.begin() + (&simplex - &simplices[0])};
        if (++counter % 10000 == 0 && simplices.size() > 100000U)
        {
          std::lock_guard<std::mutex> lock_guard(mutex);
          std::cout << "\rC++: Extracting facets ... " << counter << " / " << simplices.size();
        }
        auto first_is_face{false};
        for (auto second_it{first_it + 1}; second_it < simplices.end(); ++second_it)
        {
          first_is_face = false;
          if (simplex.is_face(*second_it))
          {
            first_is_face = true;
            break;
          }
        }
        if (!first_is_face)
        {
          std::lock_guard<std::mutex> lock_guard(mutex);
          facets.push_back(simplex.vertices());
        }
      });

  return facets;
}
