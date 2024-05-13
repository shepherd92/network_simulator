#ifndef _INFINITE_MODEL_H_
#define _INFINITE_MODEL_H_

#include "model.h"

class InfiniteNetwork;

constexpr Mark MIN_MARK{1e-2};

class InfiniteModel : virtual public Model
{
public:
    InfiniteModel(const uint32_t seed);

protected:
    float distance(const Point &first, const Point &second) const override;
};

#endif