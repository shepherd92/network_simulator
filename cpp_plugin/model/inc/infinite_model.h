#ifndef _INFINITE_MODEL_H_
#define _INFINITE_MODEL_H_

#include "model.h"

class InfiniteNetwork;

constexpr Mark MIN_MARK{1e-6};

class InfiniteModel : virtual public Model
{
public:
    std::vector<InfiniteNetwork> generate_networks(const uint32_t num_of_networks) const;

protected:
    virtual InfiniteNetwork generate_network() const = 0;
    float distance(const Point &first, const Point &second) const override;
};

#endif