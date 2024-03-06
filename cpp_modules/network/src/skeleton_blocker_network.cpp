#include "skeleton_blocker_network.h"
#include "typedefs.h"

std::vector<Dimension> SkeletonBlockerNetwork::calc_simplex_dimension_distribution()
{
    initialize_simplicial_complex_if_needed();

    ...
}

std::vector<uint32_t> SkeletonBlockerNetwork::calc_degree_sequence_simplicial_complex(
    const Dimension simplex_dimension,
    const Dimension neighbor_dimension)
{
    assert(is_valid());
    assert(neighbor_dimension > simplex_dimension);
    std::vector<uint32_t> result{};

    if (simplex_dimension == 0)
    {
        ...
    }
    else if (simplex_dimension == 1)
    {
    }
    else if (simplex_dimension == 2)
    {
    }
    else
    {
        assert(false);
    }
    std::sort(result.begin(), result.end());
}

void SkeletonBlockerNetwork::reset_simplicial_complex()
{
    skeleton_blocker_ = std::nullopt;
}

SimplexList SkeletonBlockerNetwork::convert_to_representable_simplices(const SimplexList &simplices_in) const
{
    return simplices_in;
}

bool SkeletonBlockerNetwork::is_valid() const
{
    return skeleton_blocker_.has_value();
}

void SkeletonBlockerNetwork::initialize_simplicial_complex_if_needed()
{
    if (!is_valid())
    {
        skeleton_blocker_ = SkeletonBlocker{};
    }
}