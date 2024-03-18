#ifndef _FINITE_MODEL_H_
#define _FINITE_MODEL_H_

#include "model.h"

constexpr float TORUS_SIZE{1.};
class FiniteNetwork;

class FiniteModel : virtual public Model
{
public:
    FiniteModel(const uint32_t seed);

protected:
    PointList create_points(const size_t num_of_nodes) const;
    PositionList generate_positions(const size_t num_of_points) const;
    float distance(const Point &first, const Point &second) const override;
    constexpr static float torus_size();
};

#include "finite_model.inl"

#endif