#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "degree_sequence.h"
#include "simplex.h"

std::vector<int32_t> calc_degree_sequence_interface(
    const std::vector<VertexList> &simplices_in,
    const std::vector<VertexList> &facets_in,
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    const auto simplices{select_simplices_by_dimension(create_simplices(simplices_in), simplex_dimension)};
    const auto facets{select_higher_dimensional_simplices(create_simplices(facets_in), neighbor_dimension)};

    std::vector<int32_t> degree_sequence;
    degree_sequence.reserve(simplices.size());

    auto counter{0U};
    for (const auto &simplex : simplices)
    {
        if (++counter % 100000 == 0 && simplices.size() > 1000000U)
        {
            std::cout << "\rC++: Calculating degree sequence ..."
                      << counter << " / " << simplices.size();
        }
        // Container of vertices with which the simplex forms a simplex of neighbor dimension
        SimplexSet combinations_of_remaining_vertices;

        // iterate over all facets
        for (const auto &facet : facets)
        {
            if (simplex.is_face(facet))
            {
                const auto difference{facet - simplex}; // vertices of facet that are not in the simplex
                const auto skeleton{difference.get_skeleton(neighbor_dimension - simplex_dimension)};
                std::copy(skeleton.begin(), skeleton.end(), std::inserter(combinations_of_remaining_vertices, combinations_of_remaining_vertices.end()));
            }
        }
        degree_sequence.push_back(combinations_of_remaining_vertices.size());
    }
    return degree_sequence;
}
