#ifndef _MODEL_PARAMETERS_H_
#define _MODEL_PARAMETERS_H_

#include <pybind11/numpy.h>

namespace py = pybind11;

struct ModelParameters
{
    ModelParameters(const py::array_t<double> &model_parameters_input);
    double alpha;
    double beta;
    double gamma;
    double torus_dimension;
    double torus_size;
};

#endif