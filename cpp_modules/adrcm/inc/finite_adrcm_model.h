#ifndef _FINITE_ADRCM_MODEL_H_
#define _FINITE_ADRCM_MODEL_H_

#include "adrcm_model.h"
#include "point.h"

constexpr float TORUS_SIZE{1.};

class FiniteAdrcmModel : public AdrcmModel
{
private:
    float distance(const Point &first, const Point &second) const override;
};

float FiniteAdrcmModel::distance(const Point &first, const Point &second) const
{
    const auto distance_inside{fabs(first.position() - second.position())};

    const auto distance{
        distance_inside < 0.5 * TORUS_SIZE ? distance_inside : TORUS_SIZE - distance_inside};

    return distance;
}

#endif