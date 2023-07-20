#ifndef _MODEL_PARAMETERS_H_
#define _MODEL_PARAMETERS_H_

#include <pybind11/numpy.h>

namespace py = pybind11;

struct ModelParameters
{
    ModelParameters(const py::array_t<double> &model_parameters_input);
    uint32_t max_dimension;
    uint32_t network_size;
    uint32_t num_of_nodes;
    double alpha;
    double beta;
    double gamma;
    uint32_t torus_dimension;
    double torus_size_in_1_dimension;
};

#endif