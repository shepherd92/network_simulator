#ifndef _INFINITE_ADRCM_MODEL_H_
#define _INFINITE_ADRCM_MODEL_H_

#include "adrcm_model.h"
#include "infinite_model.h"
#include "typedefs.h"

class InfiniteAdrcmModel : public InfiniteModel, public AdrcmModel
{
public:
    InfiniteAdrcmModel(const std::vector<double> &parameters_in, const uint32_t seed);
    std::vector<InfiniteNetwork> generate_networks(const uint32_t num_of_infinite_networks) const;

protected:
    InfiniteNetwork generate_network() const;

private:
    PointList create_vertices() const;
};

#endif