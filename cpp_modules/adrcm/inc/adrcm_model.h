#ifndef _ADRCM_MODEL_H_
#define _ADRCM_MODEL_H_

#include "point.h"

class AdrcmModel
{
protected:
    float torus_size() const;

private:
    virtual float distance(const Point &first, const Point &second) const = 0;
    constexpr static Dimension torus_dimension_{1};
};

#endif