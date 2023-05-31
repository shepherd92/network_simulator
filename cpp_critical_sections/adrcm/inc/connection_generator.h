#ifndef _CONNECTION_GENERATOR_H_
#define _CONNECTION_GENERATOR_H_

#include <pybind11/numpy.h>

#include "connection_generator.inl"

namespace py = pybind11;

py::array_t<int32_t> generate_finite_network_connections_default_interface(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input,
    const py::array_t<double> &model_parameters);

std::vector<py::array_t<int32_t>> generate_infinite_network_connections_default_interface(
    const py::array_t<double> &model_parameters,
    const uint32_t num_of_infinite_networks,
    const uint32_t seed);

#include "connection_generator.inl"

#endif