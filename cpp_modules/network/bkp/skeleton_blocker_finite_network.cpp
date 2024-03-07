#include <mutex>

#include "skeleton_blocker_finite_network.h"
#include "typedefs.h"

std::vector<uint32_t> SkeletonBlockerFiniteNetwork::calc_degree_sequence(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    if (!is_valid())
    {
        return Network::calc_degree_sequence(simplex_dimension, neighbor_dimension);
    }
    assert(neighbor_dimension > simplex_dimension);

    std::vector<uint32_t> result{};
    std::mutex mutex{};
    std::atomic<uint32_t> counter{0U};

    switch (simplex_dimension)
    {
    case 0:
        std::for_each(
            std::execution::seq,
            skeleton_blocker_->vertex_range().begin(),
            skeleton_blocker_->vertex_range().end(),
            [&](const VertexHandle vertex)
            {
                if (neighbor_dimension == 1)
                {
                    result.push_back(skeleton_blocker_->degree(vertex));
                }
                else
                {
                    std::atomic<uint32_t> degree{0U};
                    std::for_each(
                        std::execution::seq,
                        skeleton_blocker_->star_simplex_range(vertex).begin(),
                        skeleton_blocker_->star_simplex_range(vertex).end(),
                        [&](const SkeletonBlockerSimplex &simplex)
                        {
                            if (simplex.dimension() == neighbor_dimension)
                            {
                                ++degree;
                            }
                        });
                    result.push_back(degree);
                }
            });
        break;
    case 1:
        result = calc_degree_sequence(skeleton_blocker_->edge_range(), neighbor_dimension);
        break;
    case 2:
        result = calc_degree_sequence(skeleton_blocker_->triangle_range(), neighbor_dimension);
        break;
    default:
        assert(false);
        result = calc_degree_sequence(skeleton_blocker_->complex_simplex_range(), neighbor_dimension);
        break;
    }
    std::sort(result.begin(), result.end());
    return result;
}

template <typename SimplexIterator>
std::vector<uint32_t> SkeletonBlockerFiniteNetwork::calc_degree_sequence(
    const SimplexIterator &simplices,
    const Dimension neighbor_dimension)
{
    std::vector<uint32_t> result{};
    std::mutex mutex{};
    std::atomic<uint32_t> counter{0U};

    std::for_each(
        execution_policy,
        simplices.begin(),
        simplices.end(),
        [&](const SkeletonBlockerSimplex &simplex)
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
        });
    std::cout << "\rCalc degree sequence (skeleton blocker) ... " << counter;
    std::sort(result.begin(), result.end());
    return result;
}
