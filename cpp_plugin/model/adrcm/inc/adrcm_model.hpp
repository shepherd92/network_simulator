#ifndef _ADRCM_MODEL_HPP_
#define _ADRCM_MODEL_HPP_

#include "model.hpp"
#include "typedefs.hpp"

class AdrcmModel : virtual public Model
{
public:
    struct Parameters
    {
        Parameters(const std::vector<double> &parameters_in);
        Dimension max_dimension;
        double network_size;
        double alpha;
        double beta;
        double gamma;
    };

    AdrcmModel(const std::vector<double> &parameters_in);

protected:
    ConnectionList generate_connections(const PointList &vertices) const;

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