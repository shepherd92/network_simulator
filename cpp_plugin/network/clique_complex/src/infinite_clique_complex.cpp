#include <cassert>
#include <mutex>

#include "infinite_clique_complex.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

InfiniteCliqueComplex::InfiniteCliqueComplex(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const Mark typical_vertex_mark_,
    const MarkList &marks)
    : Network{max_dimension, vertices},
      InfiniteNetwork{typical_vertex_mark_, marks}
{
}

SimplexList InfiniteCliqueComplex::calc_neighbors(const Dimension dimension)
{
    if (dimension == 0U)
    {
        return SimplexList{};
    }
    // typical vertex is implicitly included in the neighbors
    // implement all combinations of vertices of length dimension
    Simplex simplex_all_vertices{get_vertices()};
    const auto neighbors{simplex_all_vertices.faces(dimension)};
    return neighbors;
}
