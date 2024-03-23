#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <execution>
#include <iostream>
#include <unordered_set>
#include <vector>

#ifdef DEBUG
constexpr auto execution_policy{std::execution::seq};
#else
constexpr auto execution_policy{std::execution::par_unseq};
#endif

using Dimension = int32_t;

class Point;
class Rectangle;
class Simplex;
class SimplexHash;

using PointList = std::vector<Point>;
using PointId = int32_t;
using PointIdList = std::vector<PointId>;
using InteractionId = int32_t;
using InteractionIdList = std::vector<InteractionId>;
using Mark = float;
using MarkList = std::vector<Mark>;
using Position = float;
using PositionList = std::vector<Position>;
using MarkPosition = std::pair<Mark, Position>;
using MarkPositionList = std::vector<MarkPosition>;

using ISimplexList = std::vector<PointIdList>;
using SimplexList = std::vector<Simplex>;
using SimplexSet = std::unordered_set<Simplex, SimplexHash>;
using RectangleList = std::vector<Rectangle>;
using Connection = std::pair<PointId, PointId>;
using ConnectionList = std::vector<Connection>;

#endif