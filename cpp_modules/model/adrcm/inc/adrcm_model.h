#ifndef _ADRCM_MODEL_H_
#define _ADRCM_MODEL_H_

#include <random>

#include "model.h"
#include "typedefs.h"

class AdrcmModel : virtual public Model
{
public:
    struct Parameters
    {
        Parameters(const py::array_t<double> &parameters_in);
        Dimension max_dimension;
        double network_size;
        double alpha;
        double beta;
        double gamma;
    };

    AdrcmModel(const py::array_t<double> &parameters_in);

protected:
    ConnectionList generate_connections(const PointList &vertices) const;
    SimplexList create_simplices_from_connections(const ConnectionList &connections) const;

    Dimension max_dimension() const;
    float lambda() const;
    float alpha() const;
    float beta() const;
    float gamma() const;
    bool is_default_alpha() const;

private:
    Parameters parameters_;
};

#endif