#ifndef _CONNECTION_GENERATOR_H_
#define _CONNECTION_GENERATOR_H_

#include <pybind11/numpy.h>

#include "connection_generator.inl"

namespace py = pybind11;

py::array_t<int32_t> generate_connections_default(
    const py::array_t<double> &birth_times_input,
    const py::array_t<double> &positions_input,
    const py::array_t<double> &model_parameters);

#include "connection_generator.inl"

#endif