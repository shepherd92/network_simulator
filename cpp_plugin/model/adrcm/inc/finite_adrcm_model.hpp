#ifndef _FINITE_ADRCM_MODEL_HPP_
#define _FINITE_ADRCM_MODEL_HPP_

#include "adrcm_model.hpp"
#include "finite_model.hpp"
#include "typedefs.hpp"

class FiniteCliqueComplex;

class FiniteAdrcmModel : public FiniteModel, public AdrcmModel
{
public:
    FiniteAdrcmModel(const std::vector<double> &parameters_in, const uint32_t seed);
    std::tuple<FiniteCliqueComplex, MarkPositionList> generate_network() const;
};

#endif