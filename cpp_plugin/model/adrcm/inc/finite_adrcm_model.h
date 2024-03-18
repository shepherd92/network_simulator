#ifndef _FINITE_ADRCM_MODEL_H_
#define _FINITE_ADRCM_MODEL_H_

#include "adrcm_model.h"
#include "finite_model.h"
#include "typedefs.h"

class FiniteAdrcmModel : public FiniteModel, public AdrcmModel
{
public:
    FiniteAdrcmModel(const std::vector<double> &parameters_in, const uint32_t seed);
    std::tuple<FiniteNetwork, MarkPositionList> generate_network() const;
};

#endif