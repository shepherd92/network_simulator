#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <set>
#include <vector>

#include <gudhi/Simplex_tree.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_cpp_conversion.h"

namespace py = pybind11;

class Simplex;
class Hash;

using SimplexTree = Gudhi::Simplex_tree<>;
using Dimension = uint32_t;
using VertexId = SimplexTree::Vertex_handle;
using VertexList = std::vector<VertexId>;
using ISimplex = std::vector<VertexId>;
using ISimplexList = std::vector<ISimplex>;
using SimplexList = std::vector<Simplex>;
using SimplexSet = std::unordered_set<Simplex, Hash>;

#endif