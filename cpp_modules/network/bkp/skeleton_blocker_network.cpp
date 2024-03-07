#include <mutex>

#include "skeleton_blocker_network.h"
#include "typedefs.h"

void SkeletonBlockerNetwork::add_simplex(const VertexList &simplex)
{
    initialize_simplicial_complex_if_needed();
    std::vector<VertexHandle> vertex_handles{};
    vertex_handles.reserve(simplex.size());

    std::transform(
        execution_policy,
        simplex.begin(),
        simplex.end(),
        std::back_inserter(vertex_handles),
        [&](const auto &vertex)
        {
            return VertexHandle{vertex};
        });

    std::initializer_list<VertexHandle> vertex_handles_il(vertex_handles);
    skeleton_blocker_->add_simplex(vertex_handles_il);
    SkeletonBlocker::make_complex_from_top_faces
}

SkeletonBlockerNetwork::SkeletonBlockerSimplex SkeletonBlockerNetwork::create_skeleton_blocker_simplex(
    const VertexList &vertices) const
{
    return SkeletonBlockerSimplex{};
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

void SkeletonBlockerNetwork::filter_simplicial_complex(const VertexList &vertices)
{
    assert(is_valid());
    std::for_each(
        std::execution::seq,
        vertices.begin(),
        vertices.end(),
        [&](const auto &vertex)
        {
            skeleton_blocker_->remove_vertex(VertexHandle{vertex});
        });
}