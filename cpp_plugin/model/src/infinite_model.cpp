#include "infinite_model.h"
#include "infinite_network.h"
#include "point.h"
#include "tools.h"

InfiniteModel::InfiniteModel(const uint32_t seed) : Model{seed}
{
}

float InfiniteModel::distance(const Point &first, const Point &second) const
{
    return fabs(first.position() - second.position());
}