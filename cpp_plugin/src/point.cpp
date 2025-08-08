#include "point.hpp"

MarkPositionList convert_to_mark_position_pairs(const PointList &points)
{
    MarkPositionList mark_position_pairs{};
    mark_position_pairs.reserve(points.size());
    std::transform(
        points.begin(),
        points.end(),
        std::back_inserter(mark_position_pairs),
        [](const auto &point)
        {
            return std::make_pair(point.mark(), point.position());
        });
    return mark_position_pairs;
}

PointIdList convert_to_id_list(const PointList &points)
{
    PointIdList ids{};
    ids.reserve(points.size());
    std::transform(
        points.begin(),
        points.end(),
        std::back_inserter(ids),
        [](const auto &point)
        {
            return point.id();
        });
    return ids;
}

MarkList convert_to_mark_list(const PointList &points)
{
    MarkList mark_list{};
    mark_list.reserve(points.size());
    std::transform(
        points.begin(),
        points.end(),
        std::back_inserter(mark_list),
        [](const auto &point)
        {
            return point.mark();
        });
    return mark_list;
}

void transform_points(PointList &points, const std::function<void(Point &)> &lambda)
{
    std::for_each(
        execution_policy,
        points.begin(),
        points.end(),
        [&](auto &point)
        {
            lambda(point);
        });
}