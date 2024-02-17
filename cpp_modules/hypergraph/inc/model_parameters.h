#ifndef _MODEL_PARAMETERS_H_
#define _MODEL_PARAMETERS_H_

#include <pybind11/numpy.h>

namespace py = pybind11;

struct ModelParameters
{
    ModelParameters(const py::array_t<double> &model_parameters_input);
    size_t num_of_nodes;
    size_t num_of_interactions;
    double interaction_intensity;
    double beta;
    double gamma;
    double gamma_prime;
    double torus_size;
};

#endif