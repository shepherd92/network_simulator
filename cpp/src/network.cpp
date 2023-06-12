#include "network.h"
#include "typedefs.h"

Network::Network(const simplicial_complex &simplicial_complex) : simplicial_complex_(simplicial_complex)
{
}

Complex complex(Gudhi::skeleton_blocker::make_complex_from_top_faces<Complex>(simplices.begin(), simplices.end()));

auto Network::facets() const
{
    return simplicial_complex_.maximal_simplices();
}

auto Network::simplices() const
{
    return simplicial_complex_.;
}

std::vector<int32_t> Network::calc_degree_sequence(
    const dimension simplex_dimension,
    const dimension neighbor_dimension) const
{
    const auto selected_simplices{select_simplices_by_dimension(simplex_dimension)};

    std::vector<std::vector<int32_t>> possible_neighbors;
    copy_if(facets.begin(), facets.end(),
            back_inserter(possible_neighbors),
            [neighbor_dimension](const std::vector<int32_t> &facet)
            { return facet.size() > neighbor_dimension; });

    const auto sorted_possible_neighbors{sort_simplices(possible_neighbors)};

    std::vector<int32_t> degree_sequence;

    for (const auto &simplex : selected_simplices)
    {
        // Container of vertices with which the simplex forms a simplex of neighbor dimension
        std::set<std::vector<int32_t>> combinations_of_remaining_vertices;

        // iterate over all facets
        for (const auto &facet : sorted_possible_neighbors)
        {
            // check if the facet includes the vertices of the simplex
            if (includes(facet.begin(), facet.end(), simplex.begin(), simplex.end()))
            {
                std::vector<int32_t> remaining_vertices; // vertices of facet that are not in the simplex
                std::set_difference(
                    facet.begin(), facet.end(),
                    simplex.begin(), simplex.end(),
                    std::back_inserter(remaining_vertices));
                std::sort(remaining_vertices.begin(), remaining_vertices.end());

                std::vector<int32_t> temp_combination;

                combinations(
                    remaining_vertices,
                    neighbor_dimension - simplex_dimension,
                    combinations_of_remaining_vertices,
                    temp_combination,
                    0U);
            }
        }
        degree_sequence.push_back(combinations_of_remaining_vertices.size());
    }

    return degree_sequence;
}

std::vector<simplex> Network::select_simplices_by_dimension(
    const dimension dimension)
{
    std::vector<std::vector<int32_t>> selected_simplices;

    for (const auto &simplex : simplices)
    {
        if (simplex.size() == dimension + 1)
        {
            selected_simplices.push_back(simplex);
        }
    }

    return sort_simplices(selected_simplices);
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