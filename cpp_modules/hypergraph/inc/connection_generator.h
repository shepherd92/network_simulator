#ifndef _CONNECTION_GENERATOR_H_
#define _CONNECTION_GENERATOR_H_

#include <pybind11/numpy.h>

#include "connection_generator.inl"

namespace py = pybind11;

std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>> generate_finite_network_connections_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed);

#include "connection_generator.inl"

#endif