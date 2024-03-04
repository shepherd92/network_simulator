#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <execution>
#include <iostream>
#include <set>
#include <vector>

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_cpp_conversion.h"

namespace py = pybind11;

struct SimplexTreeOptions
{
    typedef Gudhi::linear_indexing_tag Indexing_tag;
    typedef int Vertex_handle;
    typedef float Filtration_value;
    typedef uint32_t Simplex_key;
    static const bool store_key = true;
    static const bool store_filtration = false;
    static const bool contiguous_vertices = false;
};

class Simplex;
class SimplexHash;

#ifdef DEBUG
constexpr auto execution_policy{std::execution::seq};
#else
constexpr auto execution_policy{std::execution::par_unseq};
#endif

using Dimension = uint32_t;
using VertexId = int32_t;
using VertexList = std::vector<VertexId>;
using ISimplexList = std::vector<VertexList>;
using SimplexList = std::vector<Simplex>;
using SimplexSet = std::unordered_set<Simplex, SimplexHash>;

using SimplexTree = Gudhi::Simplex_tree<SimplexTreeOptions>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using PersistentCohomology = Gudhi::persistent_cohomology::Persistent_cohomology<SimplexTree, Field_Zp>;
using SimplexHandle = SimplexTree::Simplex_handle;
using SimplexHandleList = std::vector<SimplexHandle>;

#endif