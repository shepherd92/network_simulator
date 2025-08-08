#ifndef _MODEL_HPP_
#define _MODEL_HPP_

#include <random>

#include "typedefs.hpp"

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