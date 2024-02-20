#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <set>
#include <vector>

namespace py = pybind11;

class Simplex;
class Hash;

typedef uint32_t Dimension;
typedef int32_t VertexId;
typedef std::vector<VertexId> VertexList;
typedef std::vector<Simplex> SimplexList;

typedef std::unordered_set<Simplex, Hash> SimplexSet;

#endif