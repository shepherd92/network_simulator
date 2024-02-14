#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "combination.h"
#include "degree_sequence.h"
#include "simplex.h"

std::vector<int32_t> calc_degree_sequence(
    const std::vector<std::vector<int32_t>> &simplices_in,
    const std::vector<std::vector<int32_t>> &facets_in,
    const uint32_t simplex_dimension,
    const uint32_t neighbor_dimension)
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
        set_of_combinations combinations_of_remaining_vertices;

        // iterate over all facets
        for (const auto &facet : facets)
        {
            if (simplex.is_face(facet))
            {
                const auto remaining_vertices{facet - simplex}; // vertices of facet that are not in the simplex
                combination(
                    remaining_vertices,
                    neighbor_dimension - simplex_dimension,
                    combinations_of_remaining_vertices);
            }
        }
        degree_sequence.push_back(combinations_of_remaining_vertices.size());
    }
    return degree_sequence;
}
