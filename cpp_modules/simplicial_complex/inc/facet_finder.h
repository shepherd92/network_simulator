#ifndef _FACET_FINDER_
#define _FACET_FINDER_

#include <vector>

std::vector<std::vector<int32_t>> extract_facets(const std::vector<std::vector<int32_t>> &simplices);
std::vector<int32_t> calc_degree_sequence(
    const std::vector<std::vector<int32_t>> &simplices,
    const std::vector<std::vector<int32_t>> &facets,
    const uint32_t simplex_dimension,
    const uint32_t neighbor_dimension);

#endif