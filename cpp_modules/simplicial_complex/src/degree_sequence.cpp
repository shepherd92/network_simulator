#include <algorithm>
#include <iostream>
#include <execution>
#include <mutex>
#include <unordered_set>
#include <vector>

#include "degree_sequence.h"
#include "simplex.h"

std::vector<int32_t> calc_degree_sequence_interface(
    const std::vector<VertexList> &facets_in,
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    const auto facets{select_higher_dimensional_simplices(create_simplices(facets_in), neighbor_dimension)};
    std::mutex mutex;

    SimplexSet simplices{};
    std::for_each(
        std::execution::par_unseq,
        facets.begin(),
        facets.end(),
        [&](auto &&facet)
        {
            const auto skeleton{facet.get_skeleton(simplex_dimension)};
            std::lock_guard<std::mutex> lock_guard(mutex);
            for (const auto &simplex : skeleton)
            {
                simplices.insert(simplex);
            }
        });

    std::vector<int32_t> degree_sequence;
    degree_sequence.reserve(simplices.size());

    std::atomic<uint32_t> counter{0U};
    std::for_each(
        std::execution::par_unseq,
        simplices.begin(),
        simplices.end(),
        [&](auto &&simplex)
        {
            if (++counter % 100000 == 0 && simplices.size() > 1000000U)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
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
                    const auto skeleton{difference.get_skeleton(neighbor_dimension - simplex_dimension - 1)};
                    std::copy(skeleton.begin(), skeleton.end(), std::inserter(combinations_of_remaining_vertices, combinations_of_remaining_vertices.end()));
                }
            }
            std::lock_guard<std::mutex> lock_guard(mutex);
            degree_sequence.push_back(combinations_of_remaining_vertices.size());
        });

    return degree_sequence;
}
