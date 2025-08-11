#include "finite_clique_complex.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"

FiniteCliqueComplex::FiniteCliqueComplex(
    const Dimension max_dimension,
    const PointIdList &vertices,
    const ConnectionList &edges)
    : Network{max_dimension, vertices},
      FiniteNetwork{},
      CliqueComplex{edges}
{
}

FiniteCliqueComplex::FiniteCliqueComplex(FiniteCliqueComplex &&other) noexcept
    : Network{std::move(other)},
      FiniteNetwork{std::move(other)},
      CliqueComplex{std::move(other)}
{
}

// move assignment defined due to virtual base class
FiniteCliqueComplex &FiniteCliqueComplex::operator=(FiniteCliqueComplex &&other) noexcept
{
    if (this != &other)
    {
        Network::operator=(std::move(other));
        FiniteNetwork::operator=(std::move(other));
        CliqueComplex::operator=(std::move(other));
    }
    return *this;
}

SimplexList FiniteCliqueComplex::calc_simplices(const Dimension dimension)
{
    if (dimension == 0U)
    {
        std::vector<Simplex> result{};
        result.reserve(vertices_.size());
        for (const auto &vertex : vertices_)
        {
            result.emplace_back(Simplex({vertex}));
        }
        return {std::move(result)};
    }
    if (dimension == 1U)
    {
        std::vector<Simplex> result{};
        result.reserve(edges_.size());
        for (const auto &edge : edges_)
        {
            result.emplace_back(Simplex({edge.first, edge.second}));
        }
        return {std::move(result)};
    }

    assert_simplex_tree_is_built();
    SimplexList result{};

    // iterate over simplices
    for (const auto &simplex_handle : simplex_tree_->complex_simplex_range())
    {
        if (simplex_tree_->dimension(simplex_handle) == dimension)
        {
            PointIdList vertices{};
            vertices.reserve(dimension + 1U);
            for (const auto vertex : simplex_tree_->simplex_vertex_range(simplex_handle))
            {
                vertices.push_back(vertex);
            }
            result += SimplexList{std::vector{Simplex{vertices}}};
        }
    }
    return result;
}

FiniteCliqueComplex FiniteCliqueComplex::filter(const PointIdList &vertices) const
{
    return FiniteCliqueComplex{
        max_dimension_,
        vertices,
        std::move(filter_edges(vertices))};
}
