#ifndef _INFINITE_ADRCM_MODEL_H_
#define _INFINITE_ADRCM_MODEL_H_

#include "adrcm_model.h"
#include "infinite_model.h"
#include "typedefs.h"

class InfiniteAdrcmModel : public InfiniteModel, public AdrcmModel
{
public:
    InfiniteAdrcmModel(const py::array_t<double> &parameters_in, const uint32_t seed);

protected:
    InfiniteNetwork generate_network() const override;

private:
    PointList create_vertices() const;
};

#endif