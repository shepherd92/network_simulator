#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <execution>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

#ifdef DEBUG
constexpr auto execution_policy{std::execution::seq};
#else
constexpr auto execution_policy{std::execution::par_unseq};
#endif

namespace py = pybind11;

class Point;
class Rectangle;

typedef std::vector<Point> PointList;
typedef std::vector<Rectangle> RectangleList;
typedef int Id;
typedef float Mark;
typedef float Position;
typedef std::vector<Id> IdList;
typedef std::vector<Position> PositionList;
typedef std::vector<Mark> MarkList;
typedef std::pair<Id, Id> Connection;
typedef std::vector<std::pair<Id, Id>> ConnectionList;
typedef std::pair<Mark, Position> MarkPosition;
typedef std::vector<MarkPosition> MarkPositionList;
typedef std::tuple<py::array_t<Id>, py::array_t<Mark>, py::array_t<Position>> NetworkInterface;
typedef std::vector<NetworkInterface> NetworkListInterface;

#endif