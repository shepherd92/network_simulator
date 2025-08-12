#include "clique_complex.hpp"
#include "simplex.hpp"
#include "simplex_list.hpp"

CliqueComplex::CliqueComplex(const ConnectionList &edges)
    : edges_{edges}
{
}

CliqueComplex::CliqueComplex(CliqueComplex &&other) noexcept
    : edges_{std::move(other.edges_)}
{
}

// move assignment defined due to virtual base class
CliqueComplex &CliqueComplex::operator=(CliqueComplex &&other) noexcept
{
    if (this != &other)
    {
        edges_ = std::move(other.edges_);
    }
    return *this;
}

void CliqueComplex::set_edges(const ConnectionList &edges)
{
    edges_ = edges;
    reset();
}

void CliqueComplex::set_vertices(const PointIdList &vertices)
{
    Network::set_vertices(vertices);
    edges_ = std::move(filter_edges(vertices));
}

ConnectionList CliqueComplex::get_edges() const
{
    return edges_;
}

ConnectionList CliqueComplex::filter_edges(const PointIdList &vertices) const
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

void CliqueComplex::fill_simplex_tree()
{
    assert_simplex_tree_is_initialized();
    simplex_tree_->insert_batch_vertices(vertices_);
    for (const auto &edge : edges_)
    {
        // insert edges as simplices
        simplex_tree_->insert_simplex({edge.first, edge.second});
    }
    expand();
}

void CliqueComplex::expand()
{
    assert_simplex_tree_is_built();
    simplex_tree_->expansion(max_dimension_);
    reset_persistence();
}
