#include <mutex>

#include "skeleton_blocker_infinite_network.h"
#include "typedefs.h"

SkeletonBlockerInfiniteNetwork::SkeletonBlockerInfiniteNetwork(const VertexId typical_vertex_id)
    : SkeletonBlockerNetwork(),
      typical_vertex_id_{typical_vertex_id}
{
}

SkeletonBlockerInfiniteNetwork::SkeletonBlockerSimplexRange SkeletonBlockerInfiniteNetwork::get_simplices() const
{
    assert(is_valid());
    return skeleton_blocker_->star_simplex_range(typical_vertex_id_);
}

std::vector<uint32_t> SkeletonBlockerInfiniteNetwork::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    if (!is_valid())
    {
        return Network::calc_degree_sequence(simplex_dimension, neighbor_dimension);
    }
    assert(is_valid());
    assert(neighbor_dimension > simplex_dimension);

    const auto simplices{get_simplices()};

    std::vector<uint32_t> result{};
    std::mutex mutex{};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const SkeletonBlockerSimplex &simplex)
        {
            if (simplex.dimension() == simplex_dimension)
            {
                auto degree{0U};
                std::for_each(
                    std::execution::seq,
                    skeleton_blocker_->coboundary_range(simplex).begin(),
                    skeleton_blocker_->coboundary_range(simplex).end(),
                    [&](const SkeletonBlockerSimplex &coboundary_simplex)
                    {
                        if (coboundary_simplex.dimension() == neighbor_dimension)
                        {
                            ++degree;
                        }
                    });

                std::lock_guard<std::mutex> lock_guard(mutex);
                if (++counter % 1000 == 0)
                {
                    std::cout << "\rCalc degree sequence (skeleton blocker) ... " << counter;
                }
                result.push_back(degree);
            }
        });
    std::cout << "\rCalc degree sequence (skeleton blocker) ... " << counter;
    std::sort(result.begin(), result.end());
    return result;
}
