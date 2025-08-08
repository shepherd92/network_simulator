#include "finite_clique_complex.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"
#include "tools.hpp"

FiniteCliqueComplex::FiniteCliqueComplex(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ConnectionList &edges)
    : Network{max_dimension, vertices},
      FiniteNetwork{},
      edges_{edges}
{
}

SimplexList FiniteCliqueComplex::calc_simplices(const Dimension dimension)
{
    assert_simplicial_complex_is_built();
    // iterate over simplices
    SimplexList result{};

    auto &simplex_tree{simplex_tree_.value()};

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
        simplex_tree_->insert_simplex_and_subfaces({edge.first, edge.second});
    }
    expand();
}

void FiniteCliqueComplex::expand()
{
    assert_simplicial_complex_is_built();
    simplex_tree_->expansion(max_dimension_ + 1U);
    reset_persistence();
}
