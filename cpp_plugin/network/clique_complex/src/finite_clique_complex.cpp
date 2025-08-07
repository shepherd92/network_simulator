#include "finite_clique_complex.h"
#include "simplex.h"
#include "simplex_list.h"
#include "tools.h"

FiniteCliqueComplex::FiniteCliqueComplex(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ConnectionList &edges)
    : Network{max_dimension, vertices},
      FiniteNetwork{},
      edges_{edges}
{
}

FiniteCliqueComplex::FiniteCliqueComplex(const FiniteCliqueComplex &other)
    : Network{other},
      FiniteNetwork{other}
{
}

SimplexList FiniteCliqueComplex::calc_simplices(const Dimension dimension)
{
    assert_simplicial_complex_is_built();
    // iterate over simplices
    SimplexList result{};

    auto &simplex_tree{get_simplex_tree().value()};

    for (const auto &simplex_handle : simplex_tree.filtration_simplex_range())
    {
        if (simplex_tree.dimension(simplex_handle) == dimension)
        {
            PointIdList vertices{};
            vertices.reserve(dimension + 1U);
            for (const auto vertex : simplex_tree.simplex_vertex_range(simplex_handle))
            {
                vertices.emplace_back(vertex);
            }
            result += SimplexList{std::vector{Simplex{vertices}}};
        }
    }
    return result;
}

void FiniteCliqueComplex::fill_simplicial_complex()
{
    for (const auto &edge : edges_)
    {
        // insert edges as simplices
        get_simplex_tree()->insert_simplex({edge.first, edge.second});
    }
    expand();
}

void FiniteCliqueComplex::expand()
{
    assert_simplicial_complex_is_built();
    get_simplex_tree()->expansion(max_dimension_ + 1U);
    reset_persistence();
}
