#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include <pybind11/numpy.h>

namespace py = pybind11;

std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>> generate_finite_network_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed);

std::vector<std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>>> generate_infinite_networks_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t num_of_infinite_networks,
    const uint32_t seed);

#endif