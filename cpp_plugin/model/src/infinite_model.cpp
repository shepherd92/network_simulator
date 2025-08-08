#include "infinite_model.hpp"
#include "infinite_network.hpp"
#include "point.hpp"
#include "tools.hpp"

InfiniteModel::InfiniteModel(const uint32_t seed) : Model{seed}
{
}

float InfiniteModel::distance(const Point &first, const Point &second) const
{
    return fabs(first.position() - second.position());
}
