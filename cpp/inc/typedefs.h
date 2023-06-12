#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <cstdint>
#include <vector>

#include <gudhi/Skeleton_blocker.h>
#include <gudhi/Toplex_map.h>

typedef std::uint32_t vertex_id;
typedef std::uint8_t dimension;
typedef std::vector<std::pair<vertex_id, vertex_id>> connections;

typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits traits;
typedef Gudhi::skeleton_blocker::Skeleton_blocker_complex<traits> SkeletonBlocker;
typedef skeleton_blocker::Vertex_handle SkeletonBlockerVertex;
typedef skeleton_blocker::Simplex SkeletonBlockerSimplex;

typedef Gudhi::Toplex_map ToplexMap;
typedef Gudhi::Toplex_map::Simplex ToplexSimplex;

#endif