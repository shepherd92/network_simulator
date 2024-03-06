#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <execution>
#include <iostream>
#include <set>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_cpp_conversion.h"

namespace py = pybind11;

enum SimplicialComplexType : u_int8_t
{
    simplex_tree,
    skeleton_blocker
};

constexpr auto simplicial_complex_type = SimplicialComplexType::simplex_tree;

#ifdef DEBUG
constexpr auto execution_policy{std::execution::seq};
#else
constexpr auto execution_policy{std::execution::par_unseq};
#endif

using Dimension = int32_t;
using VertexId = int32_t;
using VertexList = std::vector<VertexId>;
using ISimplexList = std::vector<VertexList>;

#endif