#ifndef _MODEL_H_
#define _MODEL_H_

#include <random>

#include "typedefs.h"

class Model
{
public:
    Model(const uint32_t seed);

protected:
    MarkList generate_marks(const size_t num_nodes, const Mark min_mark = 0., const Mark max_mark = 1.) const;
    virtual float distance(const Point &first, const Point &second) const = 0;

    mutable std::mt19937 random_number_generator_;
    float determine_space_size(const PointList &points) const;

private:
};

#endif