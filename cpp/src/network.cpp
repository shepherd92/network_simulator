#include <gudhi/Skeleton_blocker.h>

#include "network.h"
#include "typedefs.h"

Network::Network(const SimplicialComplex &simplicial_complex)
    : simplicial_complex_(simplicial_complex)
{
}

Network::Network(const std::vector<Simplex> &interactions)
    : simplicial_complex_(Gudhi::skeleton_blocker::make_complex_from_top_faces<SimplicialComplex>(interactions.begin(), interactions.end()))
{
}

auto Network::facets() const
{
    for (const auto &simplex : simplices())
    {
        const auto link{simplicial_complex_.link(simplex)};
        if (link.)
        {
        }
        for (const auto &link_simplex : link.complex_simplex_range())
        {
            if (link_simplex.dimension() == codimension)
            {
                ++degree;
            }
        }

        degree_sequence.push_back(degree);
    }
}

auto Network::simplices() const
{
    return simplicial_complex_.complex_simplex_range();
}

auto Network::num_vertices() const
{
    return simplicial_complex_.num_simplices(0);
}

auto Network::num_edges() const
{
    return simplicial_complex_.num_simplices(1);
}

auto Network::num_triangles() const
{
    return simplicial_complex_.num_simplices(2);
}

auto Network::num_simplices() const
{
    return simplicial_complex_.num_simplices();
}

auto Network::calc_degree_sequence(
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

auto Network::get_simplices_by_dimension(
    const dimension dimension) const
{
    std::vector<Simplex> selected_simplices;
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