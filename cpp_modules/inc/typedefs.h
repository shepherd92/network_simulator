#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <execution>
#include <iostream>
#include <set>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_cpp_conversion.h"

namespace py = pybind11;

#ifdef DEBUG
constexpr auto execution_policy{std::execution::seq};
#else
constexpr auto execution_policy{std::execution::par_unseq};
#endif

using Dimension = int32_t;

class Point;
class Rectangle;

using PointList = std::vector<Point>;
using Mark = float;
using MarkList = std::vector<Mark>;
using Position = float;
using PositionList = std::vector<Position>;
using PointId = int32_t;
using PointIdList = std::vector<PointId>;
using ISimplexList = std::vector<PointIdList>;
using RectangleList = std::vector<Rectangle>;
using Connection = std::pair<PointId, PointId>;
using ConnectionList = std::vector<Connection>;
using MarkPosition = std::pair<Mark, Position>;
using MarkPositionList = std::vector<MarkPosition>;

typedef std::tuple<py::array_t<PointId>, py::array_t<Mark>, py::array_t<Position>> NetworkInterface;
typedef std::vector<NetworkInterface> NetworkListInterface;

#endif