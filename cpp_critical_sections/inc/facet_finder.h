#ifndef _FACET_FINDER_
#define _FACET_FINDER_

#include <vector>

using namespace std;

vector<vector<int>> extract_facets(const vector<vector<int>> &simplices);
vector<int32_t> calc_degree_sequence(
    const vector<vector<int32_t>> &simplices,
    const vector<vector<int32_t>> &facets,
    const uint32_t simplex_dimension,
    const uint32_t neighbor_dimension);

#endif