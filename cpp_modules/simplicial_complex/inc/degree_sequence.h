#ifndef _DEGREE_SEQUENCE_H_
#define _DEGREE_SEQUENCE_H_

std::vector<int32_t> calc_degree_sequence(
    const std::vector<std::vector<int32_t>> &simplices_in,
    const std::vector<std::vector<int32_t>> &facets_in,
    const uint32_t simplex_dimension,
    const uint32_t neighbor_dimension);

#endif