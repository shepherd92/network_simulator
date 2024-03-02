#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <iostream>
#include <set>
#include <vector>

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_cpp_conversion.h"

namespace py = pybind11;

using SimplexTree = Gudhi::Simplex_tree<>;
using Dimension = uint32_t;
using VertexId = SimplexTree::Vertex_handle;
using VertexList = std::vector<VertexId>;
using ISimplex = std::vector<VertexId>;
using ISimplexList = std::vector<ISimplex>;

using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using PersistentCohomology = Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree, Field_Zp>;

#endif