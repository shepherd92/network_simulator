
#include <vector>
#include <pybind11/numpy.h>

#include "numpy_cpp_conversion.h"
#include "model_parameters.h"

namespace py = pybind11;

ModelParameters::ModelParameters(const py::array_t<double> &model_parameters_input)
{
    const auto model_parameters = numpy_to_vector_1d<double>(model_parameters_input);
    max_dimension = static_cast<uint32_t>(model_parameters[0]);
    num_of_nodes = static_cast<uint32_t>(model_parameters[1]);
    alpha = model_parameters[2];
    beta = model_parameters[3];
    gamma = model_parameters[4];
    torus_dimension = static_cast<uint32_t>(model_parameters[5]);
    torus_size = model_parameters[6];
}
