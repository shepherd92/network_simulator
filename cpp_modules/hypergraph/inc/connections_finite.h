#ifndef _CONNECTIONS_FINITE_H_
#define _CONNECTIONS_FINITE_H_

#include <pybind11/numpy.h>
#include <vector>

#include "typedefs.h"

namespace py = pybind11;

class Point;

NetworkInterface generate_finite_network_interface(
    const py::array_t<double> &model_parameters_input,
    const uint32_t seed);

PointList create_points_finite(
    const size_t num_of_nodes,
    const double torus_size,
    const float exponent);

PositionList generate_positions_finite(
    const size_t num_of_nodes,
    const double torus_size);

#endif