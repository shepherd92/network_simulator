#include <execution>
#include <mutex>

#include "numpy_cpp_conversion.h"
#include "simplex.h"
#include "skeleton_extractor.h"

std::vector<py::array_t<VertexId>> create_skeleton_interface(
    const std::vector<VertexList> &simplices,
    const Dimension max_dimension)
{
    std::mutex mutex;
    std::vector<SimplexSet> simplices_by_dimension(max_dimension + 1U);

    std::for_each(
        std::execution::par_unseq,
        simplices.begin(),
        simplices.end(),
        [&](auto &&vertex_list)
        {
            const Simplex simplex{vertex_list};
            const auto skeleton{simplex.get_skeleton(max_dimension)};
            for (const auto &skeleton_simplex : skeleton)
            {
                std::lock_guard<std::mutex> lock_guard(mutex);
                simplices_by_dimension[skeleton_simplex.dimension()].insert(skeleton_simplex);
            }
        });

    std::vector<py::array_t<VertexId>> result{};
    result.reserve(max_dimension + 1U);
    for (const auto &simplices_of_current_dimension : simplices_by_dimension)
    {
        std::vector<VertexList> raw_simplex_list{};
        for (const auto &simplex : simplices_of_current_dimension)
        {
            raw_simplex_list.push_back(std::move(simplex.vertices()));
        }
        result.push_back((to_numpy(raw_simplex_list)));
    }
    return result;
}
