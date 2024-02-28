#ifndef _DEGREE_SEQUENCE_H_
#define _DEGREE_SEQUENCE_H_

#include "typedefs.h"

std::vector<int32_t> calc_degree_sequence_interface(
    const std::vector<VertexList> &facets_in,
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension);

#endif