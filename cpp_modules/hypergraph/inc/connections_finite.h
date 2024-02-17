#ifndef _CONNECTIONS_FINITE_H_
#define _CONNECTIONS_FINITE_H_

#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;

class Point;

std::tuple<py::array_t<int>, py::array_t<float>, py::array_t<float>>
generate_finite_network_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed);

std::vector<Point>
create_points_finite(
    const size_t num_of_nodes,
    const double torus_size,
    const float exponent);

std::vector<float> generate_positions_finite(
    const size_t num_of_nodes,
    const double torus_size);

#endif