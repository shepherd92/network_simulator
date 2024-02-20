#ifndef _SKELETON_EXTRACTOR_H_
#define _SKELETON_EXTRACTOR_H_

#include <vector>

#include "typedefs.h"

std::vector<py::array_t<VertexId>> create_skeleton_interface(
    const std::vector<VertexList> &simplices,
    const Dimension max_dimension);

#endif