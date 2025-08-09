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

FiniteCliqueComplex::FiniteCliqueComplex(FiniteCliqueComplex &&other) noexcept
    : Network{std::move(other)},
      FiniteNetwork{std::move(other)},
      edges_{std::move(other.edges_)}
{
}

// move assignment defined due to virtual base class
FiniteCliqueComplex &FiniteCliqueComplex::operator=(FiniteCliqueComplex &&other) noexcept
{
    if (this != &other)
    {
        Network::operator=(std::move(other));
        FiniteNetwork::operator=(std::move(other));
        edges_ = std::move(other.edges_);
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

    assert_simplicial_complex_is_built();
    SimplexList result{};

    // iterate over simplices
    for (const auto &simplex_handle : simplex_tree_->complex_simplex_range())
    {
        if (simplex_tree_->dimension(simplex_handle) == dimension)
        {
            PointIdList vertices{};
            vertices.reserve(dimension + 1U);
            std::cout << "Dimension: " << simplex_tree_->dimension(simplex_handle) << std::endl;
            std::cout << "Vertices: ";
            for (const auto vertex : simplex_tree_->simplex_vertex_range(simplex_handle))
            {
                std::cout << vertex << " ";
                vertices.push_back(vertex);
            }
            std::cout << std::endl;
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
    simplex_tree_->expansion(max_dimension_);
    reset_persistence();
}

FiniteCliqueComplex FiniteCliqueComplex::filter(const PointIdList &vertices) const
{
    return FiniteCliqueComplex{
        max_dimension_,
        vertices,
        std::move(filter_edges(vertices))};
}

void FiniteCliqueComplex::set_edges(const ConnectionList &edges)
{
    edges_ = edges;
    reset();
}

void FiniteCliqueComplex::set_vertices(const PointIdList &vertices)
{
    Network::set_vertices(vertices);
    edges_ = std::move(filter_edges(vertices));
}

ConnectionList FiniteCliqueComplex::get_edges() const
{
    return edges_;
}

ConnectionList FiniteCliqueComplex::filter_edges(const PointIdList &vertices) const
{
    ConnectionList filtered_edges;
    for (auto &edge : edges_)
    {
        if (std::find(vertices.begin(), vertices.end(), edge.first) != vertices.end() &&
            std::find(vertices.begin(), vertices.end(), edge.second) != vertices.end())
        {
            filtered_edges.push_back(edge);
        }
    }
    return filtered_edges;
}