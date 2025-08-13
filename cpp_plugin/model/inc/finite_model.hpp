#ifndef _FINITE_MODEL_HPP_
#define _FINITE_MODEL_HPP_

#include "model.hpp"

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
    static constexpr float torus_size{TORUS_SIZE};
};

#endif