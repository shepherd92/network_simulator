#include <gudhi/Skeleton_blocker.h>

#include "network.h"
#include "typedefs.h"

Network::Network(const SimplicialComplex &simplicial_complex)
    : simplicial_complex_(simplicial_complex)
{
}

Network::Network(const std::vector<Vertex> &vertices, const std::vector<Simplex> &simplices)
{
    for (const auto &vertex : vertices)
    {
        simplicial_complex_.add_simplex(Simplex{vertex});
    }
    for (const auto &simplex : simplices)
    {
        simplicial_complex_.add_simplex(simplex);
    }
}

Network::Network(const std::vector<vertex_id> &vertices, const std::vector<std::vector<vertex_id>> &interactions)
{
    for (const auto &vertex_id : vertices)
    {
        simplicial_complex_.add_simplex(Simplex{Vertex{vertex_id}});
    }
    for (const auto &interaction : interactions)
    {
        Simplex simplex{};
        for (const auto &vertex_id : interaction)
        {
            simplex.add_vertex(Vertex(vertex_id));
        }

        simplicial_complex_.add_simplex(simplex);
    }
}

void Network::add_simplices(const std::vector<std::vector<vertex_id>> &simplices)
{
    for (const auto &vertex_id_set : simplices)
    {
        // first add the vertices
        const auto max_vertex_id{*std::max_element(vertex_id_set.begin(), vertex_id_set.end())};
        if (max_vertex_id >= simplicial_complex_.num_vertices())
        {
            // some vertices are missing, add them
            const auto num_vertices_to_add{max_vertex_id - simplicial_complex_.num_vertices() + 1};
            for (auto i{0}; i < num_vertices_to_add; ++i)
            {
                simplicial_complex_.add_vertex();
            }
        }

        if (vertex_id_set.size() == 1)
        {
            continue; // all vertices were already added before
        }

        Simplex simplex{};
        for (const auto &vertex_id : vertex_id_set)
        {
            simplex.add_vertices(Vertex(vertex_id));
        }

        if (!simplicial_complex_.contains(simplex))
        {
            simplicial_complex_.add_simplex(simplex);
        }
    }
}

// std::vector<Simplex> Network::get_simplex_skeleton_for_max_dimension(const Simplex &simplex) const
// {
// }

const std::vector<Simplex> &Network::interactions()
{
    return interactions_.empty() ? facets() : interactions_;
}

void Network::interactions(const std::vector<Simplex> &interactions)
{
    interactions_ = interactions;
}

void Network::max_dimension(const dimension &max_dimension)
{
    max_dimension_ = max_dimension;
}

const std::vector<Simplex> &Network::facets()
{
    if (facets_.empty())
    {
        for (const auto &simplex : simplices())
        {
            const auto link{simplicial_complex_.link(simplex)};
            if (link.num_simplices() == 0)
            {
                facets_.push_back(simplex);
            }
        }
    }
    return facets_;
}

SimplexIterator Network::simplices() const
{
    return simplicial_complex_.complex_simplex_range();
}

uint32_t Network::num_vertices() const
{
    return simplicial_complex_.num_simplices(0);
}

uint32_t Network::num_edges() const
{
    return simplicial_complex_.num_simplices(1);
}

uint32_t Network::num_triangles() const
{
    return simplicial_complex_.num_simplices(2);
}

uint32_t Network::num_simplices() const
{
    return simplicial_complex_.num_simplices();
}

std::vector<uint32_t> Network::calc_degree_sequence(
    const dimension simplex_dimension,
    const dimension neighbor_dimension) const
{
    std::vector<uint32_t> degree_sequence{};
    const auto codimension{neighbor_dimension - simplex_dimension};
    for (const auto &simplex : simplices())
    {
        if (simplex.dimension() != simplex_dimension)
        {
            // we only care about simplices with simplex_dimension
            continue;
        }

        const auto link{simplicial_complex_.link(simplex)};
        auto degree{0U};
        for (const auto &link_simplex : link.complex_simplex_range())
        {
            if (link_simplex.dimension() == codimension)
            {
                ++degree;
            }
        }

        degree_sequence.push_back(degree);
    }

    return degree_sequence;
}

std::vector<Simplex> Network::get_simplices_by_dimension(
    const dimension dimension) const
{
    std::vector<Simplex> selected_simplices;
    if (dimension > max_dimension_)
    {
        return selected_simplices;
    }

    if (dimension == 0U)
    {
        for (const auto &vertex : simplicial_complex_.vertex_range())
        {
            selected_simplices.push_back(Simplex(vertex));
        }
    }
    else if (dimension == 1U)
    {
        for (const auto &edge : simplicial_complex_.edge_range())
        {
            selected_simplices.push_back(Simplex(Vertex(edge.m_source), Vertex(edge.m_target)));
        }
    }
    else if (dimension == 2U)
    {
        for (const auto &triangle : simplicial_complex_.triangle_range())
        {
            selected_simplices.push_back(triangle);
        }
    }
    else
    {
        for (const auto &simplex : simplices())
        {
            if (simplex.dimension() == dimension)
            {
                selected_simplices.push_back(simplex);
            }
        }
    }

    return selected_simplices;
}

void Network::combinations(
    const std::vector<int32_t> &elements,
    const uint32_t k,
    std::set<std::vector<int32_t>> &subarrays,
    std::vector<int32_t> &out,
    const uint32_t i)
{
    // do nothing for empty input
    if (elements.size() == 0)
    {
        return;
    }

    // base case: combination size is `k`
    if (k == 0)
    {
        subarrays.insert(out);
        return;
    }

    // return if no more elements are left
    if (i == elements.size())
    {
        return;
    }

    // include the current element in the current combination and recur
    out.push_back(elements[i]);
    combinations(elements, k - 1, subarrays, out, i + 1);

    // exclude the current element from the current combination
    out.pop_back(); // backtrack

    // exclude the current element from the current combination and recur
    combinations(elements, k, subarrays, out, i + 1);
}